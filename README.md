
TreeShrink is an algorithm for detecting abnormally long branches in one or more phylogenetic trees. 

- **Inputs**: 
    - One or more phylogenetic trees with branch lengths. If more than one, the trees should be on overlapping sets of species. 
    - Optional: A number `k` ≤ the total number of species.
    - Optional: a selection of one of the three implemented algorithms for outlier detection.
    - Optional: a false positive tolerance rate, α
- **Outputs**:
    - The removing sets: the set of species to be removed from each input tree to maximally reduce its diameter for each of the removal sizes 1, 2, ..., k. (TO BE UPDATED)
    - A final suggested list of species to be removed from each input tree, computed based on the selected statistical test. 
    - The shrunk trees: the input trees with the suggested leaves removed. 
    
Note that the tree diameter is the maximum distance between any two leaves of the tree. When multiple trees are available (e.g., gene trees), the statistical tests can use the information from all genes to decide what branches are too long. 

#### Publications:

An earlier version of TreeShrink is described in the following paper:

* Mai, Uyen, and Siavash Mirarab. “TreeShrink: Efficient Detection of Outlier Tree Leaves.” In RECOMB-CG 2017, Proceedings, 116–40. 2017. [doi:10.1007/978-3-319-67979-2_7](https://doi.org/10.1007/978-3-319-67979-2_7).

A journal version is currently under review. 

### Software:
The tool TreeShrink is written in Python and R. You need to have Python (either 2 or 3) and R. The tool uses the Dendropy package in Python for tree manipulation and the BMS package in R for statistical tests. TreeShrink can run on Linux, Mac OS, and Windows.

### Download:
If you have ```git```, you can simply clone the TreeShrink repository to your machine ```git clone https://github.com/uym2/TreeShrink/TreeShrink.git```. Otherwise, you can download the zip file to your machine. 

After you obtained a copy of TreeShrink, go to the TreeShrink directory. You should see the Python script ```treeshrink.py```. 

Type ```python treeshrink.py -h``` to test that TreeShrink can run on your machine.

### Installation:
All dependencies were built and included with the software deployment. If you have Python and R installed and in your PATH, no further installation is required. 

To run TreeShrink in a different location, simply add the TreeShrink directory to your PATH.

If you cannot run TreeShrink right the way, probably the included packages are incompatible with your system. Below are the clues to help you troubleshooting the problems:
1. First, please make sure that both Python and R are properly installed and are in your PATH. Type ```python``` or ```R``` to check. 
2. If you use an ```R``` version before 3.4, you probably see TreeShrink run with a warning message. Although we have not observed any problem running TreeShrink with an old ```R``` version, we recommend upgrading ```R``` to an up-to-date version. If you do not want to change your ```R``` version, we recommend rebuilding the ```BMS``` package so that it is compatible with your ```R``` version. For your convenience, we provide a script to do this.
- On a Linux/Mac OS machine, go to the TreeShrink directory and type ```bash install_BMS.sh```.
- On a Windows machine, after going to the TreeShrink directory, double click to the file ```install_BMS.cmd```. If you use command prompt, type ```install.BMS.cmd```.

### Usage: 
```bash
treeshrink.py [-h] -i INPUT [-d OUTDIR] [-t TEMPDIR] [-o OUTPUT] [-c]
                     [-k K] [-q QUANTILES] [-m MODE]
```
Arguments:
```
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input trees
  -d OUTDIR, --outdir OUTDIR
                        Output directory. Default: inferred from the input
                        trees
  -t TEMPDIR, --tempdir TEMPDIR
                        Directory to keep temporary files. If specified, the
                        temp files will be kept
  -o OUTPUT, --output OUTPUT
                        The name of the output trees. Default: inferred from
                        the input trees
  -c, --centroid        Do centroid reroot in preprocessing. Highly
                        recommended for large trees. Default: NO
  -k K, --k K           The maximum number of leaves that can be removed.
                        Default: auto-select based on the data
  -q QUANTILES, --quantiles QUANTILES
                        The quantile(s) to set threshold. Default is 0.05
  -m MODE, --mode MODE  Filtering mode: 'per-species', 'per-gene', 'all-
                        genes','auto'. Default: 'auto'
```

### Examples:
The TreeShrink package comes with several testing trees that can be found in the `test_data` folder.

The following command will produce the shrunk trees and the corresponding list of the species that were removed at false positive error rate `α = 0.05` (default)
```bash
python treeshrink.py -i test_data/mm10.trees
```

After running the command, the program will generate the folder `test_data/mm10_treeshrink/`, inside which you will find the shrunk trees (`mm10_shrunk_0.05.trees`) and the removed species (`mm10_shrunk _RS_0.05.txt`). You should see 10 trees in `mm10_shrunk_0.05.trees` corresponding to 10 trees of the input file `mm10.trees`. Accordingly, there are 10 lines in `mm10_shrunk _RS_0.05.txt`, each shows the list of species that were removed in the corresponding tree (empty lines indicating that the tree has no species removed). 

The α threshold can be adjusted using ```-q``` option. The output folder can be changed using ```-d```. Note that you can run TreeShrink with multiple α thresholds, as follow

```bash
 python treeshrink.py -i test_data/mm10.trees -q "0.05 0.10" -d test_data/mm10_treeshrink_multi
 ```
 
 The program will generate the folder `test_data/mm10_treeshrink_multi/` inside which there are two sets of shrunk trees and removing sets at α = 0.05 and α = 0.10.
 
 There are three modes in TreeShrink: 'per-gene', 'all-genes', and 'per-species'. By default TreeShrink will automatically select an appropriate mode, with highest priority to 'per-species' unless you have too few gene trees (i.e. less than 20 trees) or there are rare species (i.e. a species that occurs in less than 20 gene trees) in the dataset.
 Note that the 'auto' mode of TreeShrink never selects 'per-gene', which is only useful if the input trees are phylogenetically independent. The user has to manually select the `per-gene` mode in such a case. Use ```-m``` to change the mode.
 
```bash
python treeshrink.py -i test_data/mm10.trees -m per-gene -d test_data/mm10_treeshrink_pergene
python treeshrink.py -i test_data/mm10.trees -m all-genes -d test_data/mm10_treeshrink_allgenes
```
 
