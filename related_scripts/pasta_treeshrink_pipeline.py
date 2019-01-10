#! /usr/bin/env python
from multiprocessing.util import get_temp_dir
from sys import argv, stdout, stderr, exit
from subprocess import check_output,call
import argparse
from os.path import basename, dirname, splitext,realpath,join,normpath,isdir,isfile,exists
from os import mkdir,getcwd,rmdir,listdir
from copy import deepcopy
from shutil import rmtree, copyfile

def eprint(s):
    stderr.write(str(s)+"\n")


def main():
    
    parser = argparse.ArgumentParser()

    parser.add_argument("-i","--indir",required=True,help="The parent input directory where the trees (and alignments) can be found")
    parser.add_argument("-a","--input",required=False,help="The name of the input sequence files. Each subdirectory under it must contain a sequence file with this name. Default: input.fasta")
    parser.add_argument("-1","--ptargs1",required=False,help="A set of arguments passed to first PASTA run; e.g., '--iter-limit 1 --mask-gappy-sites=20 -d Protein'")
    parser.add_argument("-2","--tsargs",required=False,help="A set of arguments passed to tree shrink; e.g., '-q 0.01 -b 10 -s 5,2'")
    parser.add_argument("-3","--ptargs2",required=False,help="A set of arguments passed to first PASTA run; e.g., '--iter-limit 3 --mask-gappy-sites=3 -d Protein'")
    parser.add_argument("-n","--ns",required=False,help="Number of species'")
    parser.add_argument("-p","--protein",required=False,action='store_true',help="Add if inputs are proteins")
    parser.add_argument("-c","--cleanup",required=False,action='store_true',help="Cleanup the directories")
 
    args = vars(parser.parse_args())

    fn = args["input"] if args["input"] else "input.fata"
    bd = args["indir"]
    n = args["ns"] if args["ns"] else 300 
    a1 = args["ptargs1"] if args["ptargs1"] else ""
    a2 = args["tsargs"] if args["tsargs"] else ""
    a3 = args["ptargs2"] if args["ptargs2"] else ""
    
    subdirs = [d for d in listdir(bd) if exists(normpath(join(bd,d,fn)))]
    
    eprint("Found %d folders with %s" %(len(subdirs),fn))
    
    if (args["cleanup"]):
        for dir in subdirs:
            print("rm %s"%(join(bd,dir,"firstround*")))
            print("rm %s"%(join(bd,dir,"filtered*")))
        exit()
    
    #for dir in subdirs:
    
    with open('firststep.sh','w') as f:
        for dir in subdirs:
            p = ['run_pasta.py', '-i',join(bd,dir,fn),'-j' ,'firstround', a1]
            p.extend(['--iter-limit' ,'2'] if a1.find("--iter-limi") == -1 else [])
            p.extend(['--mask-gappy-sites=%d'%int(n/100)] if a1.find("--sk-gapp") == -1 else [])
            p.extend(['--num-cpus' , '1'] if a1.find("--num-cpus") == -1 else [])
            p.extend(["--max-mem-mb", "2048"]  if a1.find("--max-mem") == -1 else [])
            p.extend(["-d","Protein"] if args["protein"] else [])
            p.extend(["2>%s-error.log"%join(bd,dir,"firstround-pasta")])  
            print(" ".join(p), file=f)

    p = ["run_treeshrink.py -i", bd,"-t firstround.tre -a firstround.marker001.input.faa.aln", a2] 
    p.extend(['-b' ,'10'] if a2.find("-b") == -1 else []) 
    p.extend(['-s' ,'2,5'] if a2.find("-s") == -1 else []) 
    p.extend(['-p', "treeshrinktemps"] if a2.find("-b") == -1 else [])
    q = "0.05" if a2.find("-q") == -1 else a2.split(" ")[a2.split(" ").index("-q")+1]
    p.extend(['-q', q])
    p.extend(["|tee" ,"treeshrink-log.txt"])
    with open('secondstep.sh','w') as f:
        print(" ".join(p),file=f)
    #print("Output files written to " + outdir) 
    
    with open('thirdstep.sh','w') as f:
        for dir in subdirs:
            p = ['run_pasta.py', '-i',join(bd,dir,"firstround.marker001.input.faa_shrunk%s.aln"%q),
                 '-j' ,'filtered', 
                 '-a', '--alignment-suffix="aligned.fasta"', a3]
            p.extend(['--iter-limit' ,'1'] if a3.find("--iter-limi") == -1 else [])
            p.extend(['--num-cpus' , '1']  if a3.find("--num-cpus") == -1 else [])
            p.extend(["--max-mem-mb", "2048"]  if a3.find("--max-mem") == -1 else [])
            p.extend(["-d","Protein"] if args["protein"] else [])
            p.extend(["2>%s-error.log"%join(bd,dir,"filtered-pasta")])    
            print(" ".join(p), file=f)
    
if __name__ == "__main__":
    main()
