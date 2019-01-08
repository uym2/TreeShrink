* Version 1.3.0:
	* **Major change**: Make -s 2,5 the default. Default k was too big for many datasets. 
* Version 1.2.4:
	* Add ability to adjust default k foruma (`-s`)
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
