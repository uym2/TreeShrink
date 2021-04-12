* Version 1.3.9:
    * TreeShrink works with R 4.0
* Version 1.3.8b: (unstable)
    * Use Scipy stats instead of R to compute thresholds.
    * Experiments show that the number of species removed at each alpha value may change, but the effect of filtering on gene tree discordance remain the same.
* Version 1.3.7:
    * Fixed the error in the message when -x is specified.
* Version 1.3.6:
    * Allow exception species via -x. These species will NOT be removed in any tree.
* Version 1.3.5:
    * Read each tree instead of using Dendropy TreeList to save memory
    * Decomposition by long branches 
    * Allow mutiple-copy gene trees and gene-to-species mapping via -g
    * Output summary file
* Version 1.3.4:
    * Set recursion depth to 5000
    * Remove internal labels in pruned trees
    * Allow output prefix and allow outputing to an existing folder
    * Add --force to force overriding existing files
    * Add option to show version
    * Show help message if no argument is provided
    * Prevent conversion to scientific form
* Version 1.3.3:
	* Suppress quotes around dashes in the output trees
* Version 1.3.2:
	* Suppress quotes around underscores in the output trees
* Version 1.3.1:
	* If output directory is empty, do not complain 
* Version 1.3.0:
	* **Major change**: Make -s 2,5 the default. Default k was too big for many datasets. 
* Version 1.2.4:
	* Add ability to adjust default k formula (`-s`)
	* Output species sorted, to help comparisons
	*  New early fail checks
* Version 1.2.3:
	* Minor bug fix: make sure zero diameter doesn't cause error; just a warning. 
* Version 1.2.2:
	* Made sure .idx files are not created in the target directory. Instead, they reside in temp dir
	* Don't write alignments to temp files. 
	* Change logging messages a bit
	* Make sure all temps are in the same folder, and are all kept or removed at the end. 
* Version 1.2.1:
	* Made sure anaconda works
* Version 1.2.0:
	* **Major change**: by default, no species will be removed during the per-species test if their impact on diameter is less than 5%. Option `-b` can be used to change the default behavior.
