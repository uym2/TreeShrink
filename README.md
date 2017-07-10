# TreeShrink: efficient detection of outlier tree leaves

## Usage:

```python
treeshrink.py [-h] -i INPUT [-m METHOD] [-f FUNCTION] [-o OUTPUT] 
              [-r REMOVAL] [-q QUANTILE] [-d DIAMETER] [-g GRADIENT] [-c]
```

optional arguments:

  -h, --help            
                      
                        show this help message and exit
  
  -i INPUT, --input INPUT
  
                        input trees
                        
  -m METHOD, --method METHOD
  
                        method: ind,med,sts. Default: med
                        
  -f FUNCTION, --function FUNCTION
  
                        a function to fit to data: lnorm, kernel,lkernel. Default: lkernel
                        
  -o OUTPUT, --output OUTPUT
  
                        output trees
                        
  -r REMOVAL, --removal REMOVAL
  
                        the removing set
                        
  -q QUANTILE, --quantile QUANTILE
  
                        the cut-off quantile of the gradient to be used as threshold
                        
  -d DIAMETER, --diameter DIAMETER
  
                        list of the optimal diameter by level
                        
  -g GRADIENT, --gradient GRADIENT
  
                        list of the gradient of the diameter by level
                        
  -c, --centroid        
  
                        do centroid reroot in preprocessing
