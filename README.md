TreeShrink is an algorithm for detecting abnormally long branches in one or more phylogenetic trees. 

- **Inputs**: 
    - One or more phylogenetic trees with branch lengths. If more than one, the trees should be on overlapping sets of species (though missing data are allowed). 
    - Optional: a number `k` ≤ the total number of species.
    - Optional: a selection of one of the three implemented algorithms for outlier detection.
    - Optional: a false positive tolerance rate, α
    - Optional: a set of alignments, from which, shrunk sequences will be removed
- **Outputs**:
    - The removing list: the final suggested list of species to be removed from each input tree, computed based on the selected statistical test. 
    - The shrunk trees: the input trees with the suggested leaves removed. 
    - The filtered alignments: the input alignments (if provided) with suggested leaves removed. 
    
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
TreeShrink is written in Python and can run on Linux, Mac OS, and Windows. TreeShrink requires Python 3.8 or newer. Since v1.4.0, TreeShrink computes its statistical thresholds in pure Python and no longer requires R or the BMS R package at runtime. The DendroPy code used internally is vendored under TreeShrink's private namespace, so users do not need to install DendroPy separately. The runtime Python dependencies installed by `setup.py` are `treeswift`, `numpy`, and `scipy`.

### Anaconda
If you use Conda, try:

~~~bash
conda install treeshrink
~~~

This should work on most platforms. Let us know if it doesn't in the issues section. 

If you don't have Bioconda, you need to first add it as a channel:

~~~bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
~~~

### Install from github
If you have `git`, you can clone the TreeShrink repository to your machine `git clone https://github.com/uym2/TreeShrink.git`. Otherwise, you can download the zip file to your machine. 

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
run_treeshrink.py -h
~~~

#### FAQ

If you have trouble installing TreeShrink, below are some clues to help you troubleshoot the problems:

1. First, please make sure that Python 3.8 or newer is properly installed and is in your PATH. Type ```python``` to check.
2. Since v1.4.0, TreeShrink no longer requires R or the BMS R package. If you are using an older TreeShrink release, you may still need an R/BMS-compatible setup for threshold estimation.

## Usage: 
After installing TreeShrink, you can type 

~~~bash
run_treeshrink.py -h
~~~

to learn about all the options. 

## Examples:
The TreeShrink package comes with several testing trees that can be found in [test_data.zip](test_data.zip). If you downloaded TreeShrink from Github, you should have ```test_data.zip``` in your ```TreeShrink``` folder. If you installed using Anaconda, you should download [test_data.zip](https://github.com/uym2/TreeShrink/blob/master/test_data.zip) to your machine. Unzip ```test_data.zip``` before running the following examples.

### The simplest use case
The following command will produce the shrunk trees and the corresponding list of the species that were removed at false positive error rate `α = 0.05` (default)

~~~bash
run_treeshrink.py  -t test_data/mm10.trees
~~~

After running the command, the program will generate the folder `test_data/mm10_treeshrink/`, inside which you will find the shrunk trees (`output.trees`), the removed species (`output.txt`), and the runtime log (`output.log`). You should see 10 trees in `output.trees` corresponding to 10 trees of the input file `mm10.trees`. Accordingly, there are 10 lines in `output.txt`, each shows the list of species that were removed in the corresponding tree (empty lines indicating that the tree has no species removed). 
If you wish to customize the outputs, use ```-o``` to change the output folder and ```-O``` to change and the output prefix.

### Adjusting α threshold
The α threshold can be adjusted using ```-q``` option. 
You can run TreeShrink with multiple α thresholds, as follow

~~~bash
run_treeshrink.py  -t test_data/mm10.trees -q "0.05 0.10" -o test_data/mm10_treeshrink_multi -O shrunk
~~~
 
The program will generate the folder `test_data/mm10_treeshrink_multi/` inside which there are two sets of shrunk trees and removing sets at α = 0.05 and α = 0.10.
 
As TreeShrink is running, it prints runtime messages to the console and also writes the same messages to `<output directory>/<prefix>.log`. For the command above, the log file is `test_data/mm10_treeshrink_multi/shrunk.log`.
 
### Modes
 
 There are three modes in TreeShrink: 
 
 - 'per-gene'
 - 'all-genes'
 - 'per-species'.

By default, TreeShrink will automatically select an appropriate mode, with highest priority to 'per-species' unless you have too few gene trees (i.e. less than 20 trees) or there are rare species (i.e. a species that occurs in less than 20 gene trees) in the dataset.
 Note that the 'auto' mode of TreeShrink never selects 'per-gene', which is only useful if the input trees are phylogenetically independent. The user has to manually select the `per-gene` mode in such a case. Use ```-m``` to change the mode.
 
```bash
run_treeshrink.py  -t test_data/mm10.trees -m per-species -o test_data/mm10_treeshrink_perspecies
run_treeshrink.py  -t test_data/mm10.trees -m per-gene -o test_data/mm10_treeshrink_pergene
run_treeshrink.py  -t test_data/mm10.trees -m all-genes -o test_data/mm10_treeshrink_allgenes
```

#### Using `-i` to include alignments

The input can also be a set of alignments and trees. Alignments does not impact the outlier detection and are included only if one wishes to filter outliers from both trees and alignments. After TreeShrink detects outliers (based solely on trees), it  produces a new alignment and a new tree by removing the corresponding sequences.

To provide alignments and trees, your data must have the following structure: 

- You need to have a top folder (e.g., `test_data/mm_indir/`). 
- Inside that folder, you need to put one directory per input gene tree/alignment. In our example, these are `gene1`, `gene2`,...,`gene10`. 
- Inside each of these gene folders, you have a file for the gene tree and (optionally) a file for the alignment. All the tree files must have the same exact name. The same for all the alignment files. In our example, the trees are named `input.tree` and the alignments are named `input.fasta`.

Then, you can run TreeShrink and simply give it the name of the top folder using `-i`, the name of the gene tree files using `-t`, and the name of the alignment files (if present) using `-a`.


In this example, you will execute:

~~~ bash
run_treeshrink.py -i test_data/mm_indir -t input.tree -a input.fasta
~~~

This will produce a removing set `output.txt`, a filtered alignment `output.fasta`, and a filtered tree `output.tree` for each subdirectory `gene1`, `gene2`, ..., `gene10` of `test_data/mm_indir`. You can change the output directory using `-o` and the output prefix using `-O`.

Example input:

~~~bash
$ ls test_data/mm_indir/gene*
test_data/mm_indir/gene1:
input.fasta	input.tree

test_data/mm_indir/gene10:
input.fasta	input.tree

test_data/mm_indir/gene2:
input.fasta	input.tree

test_data/mm_indir/gene3:
input.fasta	input.tree

test_data/mm_indir/gene4:
input.fasta	input.tree

test_data/mm_indir/gene5:
input.fasta	input.tree

test_data/mm_indir/gene6:
input.fasta	input.tree

test_data/mm_indir/gene7:
input.fasta	input.tree

test_data/mm_indir/gene8:
input.fasta	input.tree

test_data/mm_indir/gene9:
input.fasta	input.tree

~~~

Example after running TreeShrink:

~~~ bash
$ run_treeshrink.py -i test_data/mm_indir -t input.tree -a input.fasta
$ ls test_data/mm_indir/gene*
test_data/mm_indir/gene1:
input.fasta	input.tree	output.fasta	output.tree	output.txt

test_data/mm_indir/gene10:
input.fasta	input.tree	output.fasta	output.tree	output.txt

test_data/mm_indir/gene2:
input.fasta	input.tree	output.fasta	output.tree	output.txt

test_data/mm_indir/gene3:
input.fasta	input.tree	output.fasta	output.tree	output.txt

test_data/mm_indir/gene4:
input.fasta	input.tree	output.fasta	output.tree	output.txt

test_data/mm_indir/gene5:
input.fasta	input.tree	output.fasta	output.tree	output.txt

test_data/mm_indir/gene6:
input.fasta	input.tree	output.fasta	output.tree	output.txt

test_data/mm_indir/gene7:
input.fasta	input.tree	output.fasta	output.tree	output.txt

test_data/mm_indir/gene8:
input.fasta	input.tree	output.fasta	output.tree	output.txt

test_data/mm_indir/gene9:
input.fasta	input.tree	output.fasta	output.tree	output.txt

~~~

## Gene folder utility

TreeShrink also includes a helper script, [`make_gene_folder.py`](make_gene_folder.py), for converting a flat directory of tree and alignment files into the per-gene folder layout expected by commands that use `-i`.

### What `make_gene_folder.py` does

The script scans one input directory, groups files by basename, and creates one subdirectory per gene or locus in a new output directory.

For each matching pair of files, it copies:

- the tree file to `input.tree`
- the alignment file to `input.fasta`

This is mainly a convenience utility for preparing data for workflows that expect this structure:

~~~bash
mydata/
  gene1/
    input.tree
    input.fasta
  gene2/
    input.tree
    input.fasta
~~~

### Basic usage

To see the command-line options:

~~~bash
python make_gene_folder.py -h
~~~

The main options are:

- `-i`, `--indir`: input directory containing tree and alignment files
- `-o`, `--outdir`: output directory to create; it must not already exist
- `-t`, `--treeExt`: tree file extension to match
- `-a`, `--alnExt`: alignment file extension to match

### Example

Suppose you start with a flat directory like this:

~~~bash
flat_data/
  gene1.tree
  gene1.fasta
  gene2.tree
  gene2.fasta
~~~

You can convert it to the folder-based layout with:

~~~bash
python make_gene_folder.py -i flat_data -o mydata -t .tree -a .fasta
~~~

This produces:

~~~bash
mydata/
  gene1/
    input.tree
    input.fasta
  gene2/
    input.tree
    input.fasta
~~~

### Notes

- The output directory must not already exist.
- Files are grouped by basename, so `gene1.tree` and `gene1.fasta` will both be placed under `mydata/gene1/`.
- Extension matching is literal against filename suffixes returned by Python's `splitext`, so values such as `.tree` and `.fasta` are the safest choices.

## Tree decomposition utility

TreeShrink also includes a helper script, [`decompose.py`](decompose.py), for splitting large phylogenetic trees into smaller subtrees. This can be useful when you want to partition a dataset into more manageable pieces before downstream analysis.

### What `decompose.py` does

The script repeatedly cuts a tree at the longest eligible branch, subject to two constraints:

- The branch length must be at least `--minBranch`.
- Both sides of the cut must contain at least `--minSize` taxa.

This process continues until no more valid cuts can be made. The output is a collection of subtrees. If an alignment is provided, the script also writes the corresponding sub-alignment for each subtree by keeping only the taxa present in that subtree.

### Basic usage

To see the command-line options:

~~~bash
python decompose.py -h
~~~

The main options are:

- `-t`, `--tree`: input tree file name or path. Default: `input.tree`
- `-i`, `--indir`: top-level input directory containing one subdirectory per gene/tree
- `-a`, `--alignment`: alignment file name present in each input subdirectory
- `--minSize`: minimum number of taxa allowed in each output subtree. Default: `20`
- `--minBranch`: minimum branch length that may be cut. Default: `1.0`
- `-o`, `--outdir`: output directory

### Single-file mode

If you provide `-t` without `-i`, the script reads trees from a single file and treats each line as one Newick tree:

~~~bash
python decompose.py -t test_data/mm10.trees --minSize 20 --minBranch 1.0
~~~

This creates an output directory named after the input tree file, for example:

~~~bash
mm10_decomposed/
~~~

Each input tree is assigned a generated name such as `gene_0001`, `gene_0002`, and so on. Every decomposed subtree is written to its own subdirectory:

~~~bash
mm10_decomposed/
  gene_0001_decomposed_1/
    tree.tre
  gene_0001_decomposed_2/
    tree.tre
  gene_0002_decomposed_1/
    tree.tre
~~~

### Directory mode with alignments

You can also decompose a collection of gene trees stored in subdirectories, optionally along with matching alignments. The expected layout is:

~~~bash
mydata/
  gene1/
    input.tree
    input.fasta
  gene2/
    input.tree
    input.fasta
~~~

Run:

~~~bash
python decompose.py -i mydata -t input.tree -a input.fasta --minSize 20 --minBranch 1.0
~~~

By default this writes results to:

~~~bash
mydata_decomposed/
~~~

For each output subtree, the script creates one directory containing:

- `tree.tre`: the decomposed subtree in Newick format
- `aln.fasta`: the alignment restricted to the taxa found in that subtree

An example output layout is:

~~~bash
mydata_decomposed/
  gene1_decomposed_1/
    tree.tre
    aln.fasta
  gene1_decomposed_2/
    tree.tre
    aln.fasta
  gene2_decomposed_1/
    tree.tre
    aln.fasta
~~~

### Notes

- The output directory must not already exist, because the script creates it with `mkdir`.
- In single-file mode, the input file should contain one Newick tree per line.
- Alignment files are optional. If `-a` is omitted, only subtree files are written.
