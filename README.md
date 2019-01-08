TreeShrink is an algorithm for detecting abnormally long branches in one or more phylogenetic trees. 

- **Inputs**: 
    - One or more phylogenetic trees with branch lengths. If more than one, the trees should be on overlapping sets of species (though missing data are allowed). 
    - Optional: A number `k` ≤ the total number of species.
    - Optional: a selection of one of the three implemented algorithms for outlier detection.
    - Optional: a false positive tolerance rate, α
    - Optiona: A set of alignments, from which, shrunk sequences will be removed
- **Outputs**:
    - The removing sets: the set of species to be removed from each input tree to maximally reduce its diameter for each of the removal sizes 1, 2, ..., k.
    - A final suggested list of species to be removed from each input tree, computed based on the selected statistical test. 
    - The shrunk trees: the input trees with the suggested leaves removed. 
    - If alignments provided, the filtered alignments with suggested leaves removed. 
    
Note that the tree diameter is the maximum distance between any two leaves of the tree. When multiple trees are available (e.g., gene trees), the statistical tests can use the information from all genes to decide what branches are too long. 

#### Publications:

The latest version of TreeShrink is described in:

* Mai, Uyen, and Siavash Mirarab. 2018. “TreeShrink: Fast and Accurate Detection of Outlier Long Branches in Collections of Phylogenetic Trees.” BMC Genomics 19 (S5): 272. https://doi.org/10.1186/s12864-018-4620-2.

An earlier version of TreeShrink is described in the following paper:

* Mai, Uyen, and Siavash Mirarab. “TreeShrink: Efficient Detection of Outlier Tree Leaves.” In RECOMB-CG 2017, Proceedings, 116–40. 2017. [doi:10.1007/978-3-319-67979-2_7](https://doi.org/10.1007/978-3-319-67979-2_7).


**Note:** Since publications, we have made two *substantial changes to defaults*,  which do impact the results:

1. We added a new option `-b` (since v1.2.0), in the per-species mode. Previously, in the per-species mode, we could get into situations where a species was removed from genes even when it was not on particularly long branches. This in because we look for outliers in the per-species distribution, and a species that usually has no impact on diameter, if it occasionally has even a small impact on the diameter, it could look like an outlier gene. This would cause removing things that shouldn't be removed. In the present version, we added an option `-b` set by default to 5 (for 5%). With this option, we can specify that a species should be removed from a  gene only if it increases the diameter by some percentage (once species with more impact are removed). We set `-b` by default to 5% currently to avoid diverging from the original publication by much. However, a higher value, like `-b 20` may make more sense for your dataset if you want to be more conservative. We suggest exploring this option. 
2. Our simple heuristic for setting `k` by default based on `n` was to use: `min(n/4,5×sqrt(n))`; this was, we now feel, too big. So, since v1.3.0, we have changed to using:  `min(n/5,2×sqrt(n))`, a value that can be adjusted using the new option `-s`.

## Installation:


#### Prerequisites:
The tool TreeShrink is written in Python and R. You need to have the following installed:

- Python (either 2 or 3) and 
- R installed in your machine. 

The tool uses the Dendropy package in Python for tree manipulation and the BMS package in R for statistical tests. TreeShrink can run on Linux, Mac OS, and Windows.


All dependencies were built and included with the software. If you have Python and R installed and in your `PATH`, no further installation is required. 

### Anaconda
If you use anaconda, try:

~~~bash
conda install -c smirarab treeshrink
~~~

this should work in most platforms. Let us know if it doesn't in the issues section. 

### Install from github
If you have `git`, you can simply clone the TreeShrink repository to your machine `git clone https://github.com/uym2/TreeShrink.git`. Otherwise, you can download the zip file to your machine. 


After downloading TreeShrink, to install, run:

~~~bash
python setup.py install
~~~
if you are not root, run:

~~~bash
python setup.py install --user
~~~

to test, run:

~~~bash
run_treeshrink -h
~~~

#### FAQ

If you cannot run TreeShrink right the way, probably the included packages are incompatible with your system. Below are some clues to help you troubleshoot the problems:

1. First, please make sure that both Python and R are properly installed and are in your PATH. Type ```python``` or ```R``` to check. 
2. If you use an ```R``` version before 3.4.0, you probably see TreeShrink run with a warning message **"package ‘BMS’ was built under R version 3.4.0"**. Although we have not observed any problem with this warning, we recommend upgrading ```R``` to version 3.4.0 or later. Alternatively, you can rebuild the ```BMS``` package so that it is compatible with your ```R``` version. For your convenience, we provide a script to do this.
	- On Linux/Mac OS machines, go to the TreeShrink directory and type ```bash install_BMS.sh```.
	- On Windows machine, after going to the TreeShrink directory, double click the file ```install_BMS.cmd```. If you use command prompt, type ```install_BMS.cmd```.



## Usage: 

```bash
run_treeshrink.py [-h] [-i INDIR] [-t TREE] [-a ALIGNMENT] [-c] [-k K]
                         [-q QUANTILES] [-m MODE] [-o OUTDIR] [-p TEMPDIR]
                         [-r LIBDIR]
```

Arguments include:

```bash
  -h, --help            show this help message and exit
  -i INDIR, --indir INDIR
                        The parent input directory where the trees (and
                        alignments) can be found
  -t TREE, --tree TREE  The name of the input tree/trees. If the input
                        directory is specified (see -i option), each
                        subdirectory under it must contain a tree with this
                        name. Otherwise, all the trees can be included in this
                        one file. Default: input.tre
  -a ALIGNMENT, --alignment ALIGNMENT
                        The name of the input alignment; can only be used when
                        the input directory is specified (see -i option). Each
                        subdirectory under it must contain an alignment with
                        this name. Default: input.fasta
  -c, --centroid        Do centroid reroot in preprocessing. Highly
                        recommended for large trees. Default: NO
  -k K, --k K           The maximum number of leaves that can be removed.
                        Default: auto-select based on the data
  -q QUANTILES, --quantiles QUANTILES
                        The quantile(s) to set threshold. Default is 0.05
  -m MODE, --mode MODE  Filtering mode: 'per-species', 'per-gene', 'all-
                        genes','auto'. Default: auto
  -o OUTDIR, --outdir OUTDIR
                        Output directory. Default: the same as input directory
                        (if it is specified) or the same as the input trees
  -p TEMPDIR, --tempdir TEMPDIR
                        Directory to keep temporary files. If specified, the
                        temp files will be kept
  -r LIBDIR, --libdir LIBDIR
                        Directory of the R libraries and scripts. Default: 3
                        layers above the current directory
```

### Examples:
The TreeShrink package comes with several testing trees that can be found in the `test_data` folder.

The following command will produce the shrunk trees and the corresponding list of the species that were removed at false positive error rate `α = 0.05` (default)

~~~bash
run_treeshrink.py  -i test_data/mm10.trees
~~~

After running the command, the program will generate the folder `test_data/mm10_treeshrink/`, inside which you will find the shrunk trees (`mm10_shrunk_0.05.trees`) and the removed species (`mm10_shrunk _RS_0.05.txt`). You should see 10 trees in `mm10_shrunk_0.05.trees` corresponding to 10 trees of the input file `mm10.trees`. Accordingly, there are 10 lines in `mm10_shrunk _RS_0.05.txt`, each shows the list of species that were removed in the corresponding tree (empty lines indicating that the tree has no species removed). 

The α threshold can be adjusted using ```-q``` option. The output folder can be changed using ```-d```. Note that you can run TreeShrink with multiple α thresholds, as follow

~~~bash
run_treeshrink.py  -i test_data/mm10.trees -q "0.05 0.10" -d test_data/mm10_treeshrink_multi
~~~
 
 The program will generate the folder `test_data/mm10_treeshrink_multi/` inside which there are two sets of shrunk trees and removing sets at α = 0.05 and α = 0.10.
 
As TreeShrink is running, it will output a bunch of messages to the console. We suggest saving these in a text log file:

~~~bash
run_treeshrink.py  -i test_data/mm10.trees > test_data/mm10.trees.treeshrinklog.txt
~~~
 
### Modes
 
 There are three modes in TreeShrink: 
 
 - 'per-gene'
 - 'all-genes'
 - 'per-species'.

By default, TreeShrink will automatically select an appropriate mode, with highest priority to 'per-species' unless you have too few gene trees (i.e. less than 20 trees) or there are rare species (i.e. a species that occurs in less than 20 gene trees) in the dataset.
 Note that the 'auto' mode of TreeShrink never selects 'per-gene', which is only useful if the input trees are phylogenetically independent. The user has to manually select the `per-gene` mode in such a case. Use ```-m``` to change the mode.
 
```bash
run_treeshrink.py  -i test_data/mm10.trees -m per-species -d test_data/mm10_treeshrink_perspecies
run_treeshrink.py  -i test_data/mm10.trees -m per-gene -d test_data/mm10_treeshrink_pergene
run_treeshrink.py  -i test_data/mm10.trees -m all-genes -d test_data/mm10_treeshrink_allgenes
```

#### Using `-i` to input alignments

The input can also be a set of alignments and trees. Alignments will no impact the outlier detection. The only way in which they are used is the following. After TreeShrink detects (based on trees) what sequences are outliers, it then removes those sequences from both the tree and the alignment and produces a new filtered alignment and tree.  

To provide alignment and trees, your data need to follow a certain structure, as shown in the example below. 

- You need to a top folder (e.g., `allgenes` below). 
- Inside that folder, you need to put one directory per input gene tree/alignment. In our example, these are `4048`, `4103`, `4218`, `4234`, and `4308`. 
- Inside each of these gene folders, you have a file for the gene tree and all the tree files should have the same exact name. In this case, the name is `RAxML_bipartitions.final `. 
- Inside each of these gene folders, you can have a file for the gene alignment and all the alignment files should have the same exact name. In this case, the name is `alignment.fasta`.
- Then, you can run TreeShrink and simply give it the name of the top folder using `-i`, the name of the gene tree files using `-t` and the name of the alignment files (if present) using `-a`. In this example, you will exectue:

~~~ bash
run_treeshrink.py -i allgenes -t RAxML_bipartitions.final -a alignment.fasta
~~~

- This will produce the removing sets, in addition to a filtered alignment and a filtered tree, as shown in the second list below. Note the addition of files for the shrunk tree (`RAxML_bipartitions_shrunk_0.05.final`) and the shrunk alignment (`alignment_shrunk0.05.fasta`) in addition to the files that give the name of filtered sequences (`RAxML_bipartitions_shrunk_RS_0.05.txt`).

Example input:

~~~bash
> ls allgenes/*

allgenes/4048:
RAxML_bipartitions.final
alignment.fasta

allgenes/4103:
RAxML_bipartitions.final
alignment.fasta

allgenes/4218:
RAxML_bipartitions.final
alignment.fasta

allgenes/4234:
RAxML_bipartitions.final
alignment.fasta

allgenes/4308:
RAxML_bipartitions.final
alignment.fasta

~~~


Example after running TreeShrink:

~~~ bash
allgenes/4048:
RAxML_bipartitions.final
RAxML_bipartitions_shrunk_0.05.final
RAxML_bipartitions_shrunk_RS_0.05.txt
alignment.fasta
alignment_shrunk0.05.fasta

allgenes/4103:
RAxML_bipartitions.final
RAxML_bipartitions_shrunk_0.05.final
RAxML_bipartitions_shrunk_RS_0.05.txt
alignment.fasta
alignment_shrunk0.05.fasta

allgenes/4218:
RAxML_bipartitions.final
RAxML_bipartitions_shrunk_0.05.final
RAxML_bipartitions_shrunk_RS_0.05.txt
alignment.fasta
alignment_shrunk0.05.fasta

allgenes/4234:
RAxML_bipartitions.final
RAxML_bipartitions_shrunk_0.05.final
RAxML_bipartitions_shrunk_RS_0.05.txt
alignment.fasta
alignment_shrunk0.05.fasta

allgenes/4308:
RAxML_bipartitions.final
RAxML_bipartitions_shrunk_0.05.final
RAxML_bipartitions_shrunk_RS_0.05.txt
alignment.fasta
alignment_shrunk0.05.fasta

~~~
 
