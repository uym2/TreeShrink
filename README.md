
TreeShrink is an algorithm for detecting abnormally long branches in one or more phylogenetic trees. 

- **Inputs**: 
    - One or more phylogenetic trees with branch lengths. If more than one, the trees should be on overlapping sets of species. 
    - Optional: A number `k` ≤ the total number of species.
    - Optional: a selection of one of the three implemented algorithms for outlier detection.
    - Optional: a false positive tolerance rate, α
- **Outputs**:
    - The removing sets: the set of species to be removed from each input tree to maximally reduce its diameter for each of the removal sizes 1, 2, ..., k. 
    - A final suggested list of species to be removed from each input tree, computed based on the selected statistical test. 
    - The shrunk trees: the input trees with the suggested leaves removed. 
    
Note that the tree diameter is the maximum distance between any two leaves of the tree. When multiple trees are available (e.g., gene trees), the statistical tests can use the information from all genes to decide what branches are too long. 

#### Publications:

An earlier version of TreeShrink is described in the following paper:

* Mai, Uyen, and Siavash Mirarab. “TreeShrink: Efficient Detection of Outlier Tree Leaves.” In RECOMB-CG 2017, Proceedings, 116–40. 2017. [doi:10.1007/978-3-319-67979-2_7](https://doi.org/10.1007/978-3-319-67979-2_7).

A journal version is currently under review. 

### Installation:
(TO BE UPDATED)

To be able to use TreeShrink, you need to first install the following packages. 

**Dependencies:**

- Python and [Dendropy (version 4.2.0 recommended)](https://pythonhosted.org/DendroPy/downloading.html)
- R and the [BMS package](http://bms.zeugner.eu/getBMS/)

After installing these, you need to [TO BE UPDATED]. 

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
                        genes'. Default: 'per-species'
```

### Examples:
The TreeShrink package comes with several testing trees that can be found in the `test_data` folder.

The following command will produce the shrunk trees and the corresponding removing sets at false positive error rate `α = 0.05` (default)
```bash
python treeshrink.py -i test_data/mm10.trees
```

After running the command, the program will generate the folder `test_data/mm10_treeshrink/`, inside which you will find the shrunk trees (`mm10_shrunk_0.05.trees`) and the removing sets (`mm10_shrunk _RS_0.05.txt`).

The α threshold can be adjusted using ```-q``` option. The output folder can be changed using ```-d```. Note that you can run TreeShrink with multiple α thresholds, as follow

```bash
 python treeshrink.py -i test_data/mm10.trees -q "0.05 0.10" -d test_data/mm10_treeshrink_multi
 ```
 
 The program will generate the folder `test_data/mm10_treeshrink_multi/` inside which there are two sets of shrunk trees and removing sets at α = 0.05 and α = 0.10.
 
 The default mode of TreeShrink is "per-species", which is designed to find outliers for a collection of phylogenetic trees. In this mode, the statistical tests are performed for each species. We recommend switching to the "all-genes" mode if there are rare species in the dataset. Besides, if the input trees are not phylogenetically dependent, one should use the "per-gene" mode instead. Use ```-m``` to change the mode.
 
```bash
python treeshrink.py -i test_data/mm10.trees -m per-gene
python treeshrink.py -i test_data/mm10.trees -m all-genes
```
 
