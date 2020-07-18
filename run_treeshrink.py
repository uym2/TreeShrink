#! /usr/bin/env python

from treeshrink.sequence_lib import sample_from_list
import treeshrink
from treeshrink.optimal_filter_lib import TreeFilter
from treeshrink.tree_lib import prune_tree, get_taxa,tree_as_newick
from sys import argv, stdout,setrecursionlimit
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
import re
    
def make_dir(dirName):
    if exists(dirName) and isdir(dirName):
        return False
    mkdir(dirName)
    return True

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i","--indir",required=False,help="The parent input directory where the trees (and alignments) can be found")
    parser.add_argument("-t","--tree",required=False,help="The name of the input tree/trees. If the input directory is specified (see -i option), each subdirectory under it must contain a tree with this name. Otherwise, all the trees can be included in this one file. Default: input.tre")
    parser.add_argument("-a","--alignment",required=False,help="The name of the input alignment; can only be used when the input directory is specified (see -i option). Each subdirectory under it must contain an alignment with this name. Default: input.fasta")
    parser.add_argument("-c","--centroid",required=False,action='store_true',help="Do centroid reroot in preprocessing. Highly recommended for large trees. Default: NO")
    parser.add_argument("-k","--k",required=False,help="The maximum number of leaves that can be removed. Default: auto-select based on the data; see also -s")
    parser.add_argument("-s","--kscaling",required=False,help="If -k not given, we use k=min(n/a,b*sqrt(n)) by default; using this option, you can set the a,b constants; Default: '5,2'")
    parser.add_argument("-q","--quantiles",required=False,help="The quantile(s) to set threshold. Default is 0.05")
    parser.add_argument("-b","--minImpact",required=False,help="Do not remove species on the per-species test if their impact on diameter is less than x%% where x is the given value. Default: 5")
    parser.add_argument("-m","--mode",required=False,help="Filtering mode: 'per-species', 'per-gene', 'all-genes','auto'. Default: auto")
    parser.add_argument("-o","--outdir",required=False,help="Output directory. Default: If the input directory is specified, outputs will be placed in that input directory. Otherwise, a directory with the suffix 'treeshrink' will be created in the same place as the input trees")
    parser.add_argument("-O","--outprefix",default="output",required=False,help="Output name prefix. Default: 'output'")
    parser.add_argument("-p","--tempdir",required=False,help="Directory to keep temporary files. If specified, the temp files will be kept")
    parser.add_argument("-v","--version",required=False,action='store_true',help="Show TreeShrink version.")

    if len(argv) == 1:
        parser.print_help()
        exit(0)

    args = vars(parser.parse_args())
    
    if args["version"]:
        print(treeshrink.PROGRAM_VERSION)
        exit(0)

    setrecursionlimit(5000)

    print("Launching " + treeshrink.PROGRAM_NAME + " version " + treeshrink.PROGRAM_VERSION)
    print(treeshrink.PROGRAM_NAME + " was called as follow")
    print(" ".join(argv))


    MIN_OCC = 20
    MIN_TREE_NUM = 20

    #libdir = args["libdir"] if args["libdir"] else dirname(dirname(realpath(treeshrink.__file__)))
    libdir = dirname(dirname(realpath(treeshrink.__file__)))

    tempdir = set_tmp_dir(args["tempdir"])  
    
    quantiles = [ q for q in args["quantiles"].split()] if args["quantiles"] else ["0.05"]
    
    minImpact = (float(args["minImpact"])/100)+1 if args["minImpact"] else 1.05
    
    scaling = [int(x) for x in args["kscaling"].split(",")] if  args["kscaling"] else [5,2]

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
    elif args["indir"]:
        outdir = args["indir"]
    else:
        outdir = splitext(intrees)[0] + "_treeshrink"
    if not make_dir(outdir):
        print("Warning: outputs will be written to an existing directory " + outdir)

    ''' Check to make sure output can be written
    if args["indir"]:
        i = 0
        fName,ext = splitext(basename(intrees))
        for sd in subdirs:
            outfile = normpath(join(outdir,sd, fName + "_shrunk_RS_" + quantiles[i] + ".txt"))
            with open(outfile,'w') as f:
                pass '''


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
        a_filter = TreeFilter(ddpTree=a_tree,centroid_reroot=args["centroid"],scaling=scaling)
        a_filter.optFilter(d=k)

        # compute species feature (i.e. the max ratio associated with each species for this gene tree)
        mapping = {}
        #print(a_filter.min_diams)
        for i in range(1,len(a_filter.min_diams)):
            if a_filter.min_diams[i] == 0:
                print("Warning: tree %d has no diameter (has only zero branch lengths) after removing %d sequences." %(t+1,i))
                break
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
            print("TreeShrink will run in 'Per-species' mode ...    ")

# fit kernel density to the per-species distributions and compute per-species threshold (per-species mode)
    if mode == 'per-species':
        for s in sorted(species_map):
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
                thresholds[i] = max(minImpact,float(check_output(["Rscript",normpath(join(libdir,"R_scripts","find_threshold_lkernel.R")),libdir,filename,q]).lstrip().rstrip()[5:]))
                print("%s:\n\t will be cut in %d trees where its impact is above %f for quantile %s" %(s,sum(1 for x in species_map[s] if x>thresholds[i]),thresholds[i],q,))
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

    print("Writing output ...\n")
# Dendropy's filter_leaf_nodes() seems to have problem
# i.e. it produces the trees that the treecmp tool cannot compute the MS distance (need further exploration)
# use home-made code to prune the tree instead

    fName,ext = splitext(basename(args["tree"]))
    prefix = args["outprefix"]
    counter = 0
    # check if the outdir already has files with the specified prefix
    for File in listdir(outdir):
        if File.startswith(prefix):
             search_counter = re.search(r'\d+', File[len(prefix):])
             counter = max(counter,1 if not search_counter else int(search_counter.group())+1)
    if counter > 0:
        prefix = prefix + str(counter)
     
    for i,RS in enumerate(removing_sets):
        trees_shrunk = deepcopy(trees)
        RS_tag = '' if (len(removing_sets) < 2) else '_RS_shrunk_' + quantiles[i]
        tree_tag = '' if (len(removing_sets) < 2) else '_tree_shrunk_' + quantiles[i]
        aln_tag = '' if (len(removing_sets) < 2) else '_aln_shrunk_' + quantiles[i]
        
        if args["indir"] is None:
            outfile = normpath(join(outdir,prefix + RS_tag + ".txt"))
            with open(outfile,'w') as f:
                for item in RS:
                    for s in item:
                        f.write(s + "\t")
                    f.write("\n")
            for tree,rs in zip(trees_shrunk,RS):
                prune_tree(tree,rs)
                tree_as_newick(tree,outfile=normpath(join(outdir,prefix + tree_tag + ext)),append=True)

            #trees_shrunk.write_to_path(normpath(join(outdir,fName + "_" + quantiles[i] + ext)),'newick',unquoted_underscores=True,real_value_format_specifier=".16g")  
        else:
            for sd,item in zip(subdirs,RS):
                outfile = normpath(join(outdir,sd, prefix + RS_tag + ".txt"))
                with open(outfile,'w') as f:
                    for s in item:
                        f.write(s + "\t")
            for sd,tree,rs in zip(subdirs,trees_shrunk,RS):
                L = set(x.taxon.label for x in tree.leaf_node_iter())
                prune_tree(tree,rs)
                treefile = normpath(join(outdir,sd, prefix + tree_tag + ext))
                #tree.write_to_path(treefile,'newick',unquoted_underscores=True,real_value_format_specifier=".16g")
                tree_as_newick(tree,outfile=treefile,append=False)
                
                aln_filename = args["alignment"] if args["alignment"] else "input.fasta"
                alnName,alnExt = splitext(aln_filename)
                input_aln = normpath(join(args["indir"],sd,aln_filename))
                if isfile(input_aln): 
                    output_aln = normpath(join(outdir,sd,prefix+aln_tag+alnExt))
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
