#! /usr/bin/env python
from multiprocessing.util import get_temp_dir

def main():
    from treeshrink.sequence_lib import sample_from_list
    import treeshrink
    from treeshrink.optimal_filter_lib import TreeFilter
    from treeshrink.tree_lib import prune_tree, get_taxa
    from sys import argv, stdout
    from math import sqrt
    from subprocess import check_output,call
    import argparse
    from dendropy import Tree, TreeList
    from os.path import basename, dirname, splitext,realpath,join,normpath,isdir,isfile,exists
    from os import mkdir,getcwd,rmdir,listdir
    from copy import deepcopy
    from shutil import rmtree, copyfile
    from treeshrink.alignment import CompactAlignment
    from treeshrink import set_tmp_dir, get_tmp_dir, get_tmp_file

    print("Launching " + treeshrink.PROGRAM_NAME + " version " + treeshrink.PROGRAM_VERSION)
    

    parser = argparse.ArgumentParser()

    parser.add_argument("-i","--indir",required=False,help="The parent input directory where the trees (and alignments) can be found")
    parser.add_argument("-t","--tree",required=False,help="The name of the input tree/trees. If the input directory is specified (see -i option), each subdirectory under it must contain a tree with this name. Otherwise, all the trees can be included in this one file. Default: input.tre")
    parser.add_argument("-a","--alignment",required=False,help="The name of the input alignment; can only be used when the input directory is specified (see -i option). Each subdirectory under it must contain an alignment with this name. Default: input.fasta")
    parser.add_argument("-c","--centroid",required=False,action='store_true',help="Do centroid reroot in preprocessing. Highly recommended for large trees. Default: NO")
    parser.add_argument("-k","--k",required=False,help="The maximum number of leaves that can be removed. Default: auto-select based on the data")
    parser.add_argument("-q","--quantiles",required=False,help="The quantile(s) to set threshold. Default is 0.05")
    parser.add_argument("-b","--minimpact",required=False,help="Do not remove species on the per-species test if their impact on diameter is less than MINIPACT%% where x is the given value. Default: 5")
    parser.add_argument("-m","--mode",required=False,help="Filtering mode: 'per-species', 'per-gene', 'all-genes','auto'. Default: auto")
    parser.add_argument("-o","--outdir",required=False,help="Output directory. Default: the same as input directory (if it is specified) or the same as the input trees")
    parser.add_argument("-p","--tempdir",required=False,help="Directory to keep temporary files. If specified, the temp files will be kept")
    parser.add_argument("-r","--libdir",required=False,help="Directory of the R libraries and scripts. Default: 2 layers above the treeshrink package")

    args = vars(parser.parse_args())


    MIN_OCC = 20
    MIN_TREE_NUM = 20

    libdir = args["libdir"] if args["libdir"] else dirname(dirname(realpath(treeshrink.__file__)))

    tempdir = set_tmp_dir(args["tempdir"])  
    
    quantiles = [ q for q in args["quantiles"].split()] if args["quantiles"] else ["0.05"]
    
    minimpact = (float(args["minimpact"])/100)+1 if args["minimpact"] else 1.05

    if args["indir"]:
        treename = splitext(args["tree"])[0] if args["tree"] else "input"
        subdirs = [d for d in listdir(args["indir"]) if exists(normpath(join(args["indir"],d,args["tree"] if args["tree"] else "input.tre")))]
        intrees = get_tmp_file(treename + ".trees")
        with open(intrees,'w') as fout:
            for d in subdirs:
                treename = args["tree"] if args["tree"] else "input.tre"
                treefile = normpath(join(args["indir"],d,treename))
                if exists(treefile):
                    fout.write(open(treefile,'r').read())                
    else:
        intrees = args["tree"]


    mode = args["mode"] if args["mode"] else 'auto'

    k = int(args["k"]) if args["k"] else None

    if args["outdir"]:
        outdir = args["outdir"] 
        mkdir(outdir)
    elif args["indir"]:
        outdir = args["indir"]
    else:
        outdir = splitext(intrees)[0] + "_treeshrink"
        mkdir(outdir)

    trees = TreeList.get(path=intrees,schema='newick',preserve_underscores=True)

    if mode=='auto' and len(trees) < MIN_TREE_NUM:
        print("There are only " + str(len(trees)) + " gene trees in the dataset.")
        print("TreeShrink will run in 'All-genes' mode")
        mode='all-genes'

    gene_list = [[] for i in range(len(trees))]
    species_map = {}
    occ = {}
    removing_sets = [ [ [ ] for i in range(len(trees)) ] for j in range(len(quantiles)) ]

    for t,a_tree in enumerate(trees):
        # solve k-shrink
        a_filter = TreeFilter(ddpTree=a_tree,centroid_reroot=args["centroid"])
        a_filter.optFilter(d=k)

        # compute species feature (i.e. the max ratio associated with each species for this gene tree)
        mapping = {}
        for i in range(1,len(a_filter.min_diams)):
            r = a_filter.min_diams[i-1]/a_filter.min_diams[i]
            removals = a_filter.list_removals(d=i)
            for s in removals:
                mapping[s] = r if s not in mapping else max(mapping[s],r)
        
        # gather per-species distributions and per-gene species features
        for s in mapping:
            if mode == 'per-species' or mode == 'auto':
                species_map[s] = [mapping[s]] if s not in species_map else species_map[s]+[mapping[s]]
            if mode == 'per-species' or mode == 'all-genes' or mode == 'auto':
                gene_list[t].append((s,mapping[s]))
        
        # fit kernel density to this gene's species features (per-gene mode)
        if mode == 'per-gene':
            filename = get_tmp_file("gene_%s.dat" %str(t))
            with open(filename,'w') as f:
                for s in mapping:
                    f.write(str(mapping[s]))
                    f.write("\n")
                #n_missing = len(list(a_tree.leaf_node_iter())) - len(mapping)
                #for i in range(n_missing):
                #    f.write("1.0")
                #    f.write("\n")
            if len(mapping) > 1:
                for i,q in enumerate(quantiles):
                    threshold = float(check_output(["Rscript",normpath(join(libdir,"R_scripts","find_threshold_loglnorm.R")),filename,q]).lstrip().rstrip()[4:]) 
                    #print("Threshold: ", threshold)
                    for s in mapping:
                        if mapping[s] > threshold: 
                            removing_sets[i][t].append(s)
        # update taxon occupancy (only for per-species mode)
        if mode == 'per-species' or mode == 'auto':
            for n in a_tree.leaf_node_iter():
                s = n.taxon.label
                occ[s] = 1 if not s in occ else occ[s]+1
    
    if mode == 'auto' or mode == 'per-species':
        flag = False
        for s in occ:
            if occ[s] < MIN_OCC:
                print ("Species " + s + " only exists in " + str(occ[s]) + " gene trees")
                flag = True
        if flag:
            if mode == 'auto':
                mode = 'all-genes'
                print ("There are species with low occupancy in the dataset. TreeShrink will run in 'All-genes' mode")
            else:
                print ("WARNING: 'Per-species' mode was selected for a dataset having low occupancy species. Consider switching to 'All-genes' mode")
        elif mode == 'auto':
            mode = 'per-species'
            print("Finish preprocessing. TreeShrink will run in 'Per-species' mode ...    ")

# fit kernel density to the per-species distributions and compute per-species threshold (per-species mode)
    if mode == 'per-species':
        for s in species_map:
            l = len(species_map[s])
            for i in range(occ[s]-l):
                species_map[s].append(1)
            filename = get_tmp_file(s + ".dat")
            with open(filename,'w') as f:
                for v in species_map[s]:
                    f.write(str(v))
                    f.write("\n")
            thresholds = [ 0 for i in range(len(quantiles)) ]        
            for i,q in enumerate(quantiles): 
                thresholds[i] = max(minimpact,float(check_output(["Rscript",normpath(join(libdir,"R_scripts","find_threshold_lkernel.R")),libdir,filename,q]).lstrip().rstrip()[5:]))
                print("%s:\n\t will be cut in %d tree where its impact is above %f for quantile %s" %(s,sum(1 for x in species_map[s] if x>thresholds[i]),thresholds[i],q,))
            species_map[s] = (species_map[s],thresholds)

        for t,gene in enumerate(gene_list):
            for s,r in gene:
                for i,threshold in enumerate(species_map[s][1]):
                    if r > threshold:
                        removing_sets[i][t].append(s)
                    

# fit kernel density to all the species features across all genes and compute the global threshold (all-gene mode) 
    if mode == 'all-genes':
        filename = get_tmp_file("all_genes" + ".dat")
        with open(filename,'w') as f:
            for gene in gene_list:
                for s,r in gene:
                    f.write(str(r))
                    f.write("\n")
        for i,q in enumerate(quantiles):
            threshold = float(check_output(["Rscript",normpath(join(libdir,"R_scripts","find_threshold_lkernel.R")),libdir,filename,q]).lstrip().rstrip()[5:])
            for t,gene in enumerate(gene_list):
                for s,r in gene:
                    if r > threshold:
                        removing_sets[i][t].append(s)

    print("Writing output ...\nNumber of sequences cut from each tree: ")
# Dendropy's filter_leaf_nodes() seems to have problem
# i.e. it produces the trees that the treecmp tool cannot compute the MS distance (need further exploration)
# use home-made code to prune the tree instead

    #treeName,treeExt = splitext(basename(intrees))
    #outtrees = args["output"] if args["output"] else treeName + "_shrunk" + treeExt
    fName,ext = splitext(basename(intrees))
    
    for i,RS in enumerate(removing_sets):
        trees_shrunk = deepcopy(trees)
        
        if args["indir"] is None:
            outfile = normpath(join(outdir,fName + "_RS_" + quantiles[i] + ".txt"))
            with open(outfile,'w') as f:
                for item in RS:
                    for s in item:
                        f.write(s + "\t")
                    f.write("\n")
            for tree,rs in zip(trees_shrunk,RS):
                prune_tree(tree,rs)
                print(len(rs),end="\t")
            trees_shrunk.write_to_path(normpath(join(outdir,fName + "_" + quantiles[i] + ext)),'newick')  
            print() 
        else:
            for sd,item in zip(subdirs,RS):
                outfile = normpath(join(outdir,sd, fName + "_shrunk_RS_" + quantiles[i] + ".txt"))
                print("%s: %s" %(sd,len(item)),end="\t")
                with open(outfile,'w') as f:
                    for s in item:
                        f.write(s + "\t")
            print()            
            for sd,tree,rs in zip(subdirs,trees_shrunk,RS):
                L = set(x.taxon.label for x in tree.leaf_node_iter())
                prune_tree(tree,rs)
                treeName,treeExt = splitext(args["tree"])
                treefile = normpath(join(outdir,sd, treeName + "_shrunk_" + quantiles[i] + treeExt))
                tree.write_to_path(treefile,'newick',unquoted_underscores=True,real_value_format_specifier=".16g")
                
                aln_filename = args["alignment"] if args["alignment"] else "input.fasta"
                alnName,alnExt = splitext(aln_filename)
                input_aln = normpath(join(args["indir"],sd,aln_filename))
                if isfile(input_aln): 
                    output_aln = normpath(join(outdir,sd,alnName+"_shrunk"+quantiles[i]+alnExt))
                    alg = CompactAlignment()
                    alg.read_file_object(input_aln,'fasta')
                    S=set(alg.keys())
                    if (L.difference(alg.keys())) or S.difference(L):
                        print("ERROR: For gene %s, alignment names don't match tree names. Will skip it.\n\tonly in tree:\t%s\n\tonly in alignment:\t%s"%(sd,str(L.difference(S)),str(S.difference(L))))
                    else:
                        alg.remove_all(rs)
                        alg.mask_gapy_sites(1)
                        alg.write(output_aln,'fasta')

    if not args["tempdir"]:
        rmtree(tempdir)
#    call(["rm","-r",tempdir])

    print("Output files written to " + outdir) 

    
if __name__ == "__main__":
    main()
