# Threshold Library Implementation Note

`treeshrink/threshold_lib.py` replaces the R-based threshold calculations used by earlier TreeShrink versions. TreeShrink v1.4.0 no longer requires `Rscript` or the R `BMS` package at runtime.

## Integration

`run_treeshrink.py` imports:

- `threshold_l_kernel` for all-genes and per-species thresholding
- `threshold_loglnorm` for per-gene thresholding

These replace subprocess calls to:

- `R_scripts/find_threshold_lkernel.R`
- `R_scripts/find_threshold_loglnorm.R`

## Implemented Methods

`threshold_l_kernel` mirrors the previous R kernel-density threshold:

```r
x = y[y > 0]
exp(quantile.density(density(log(x), adjust = 1), p = 1 - e))
```

The Python implementation filters positive finite values, applies `log`, uses R's `bw.nrd0` bandwidth rule, computes a Gaussian kernel density, interpolates the density quantile, exponentiates, and rounds to six decimals.

`threshold_loglnorm` mirrors the previous log-lognormal path by fitting a normal distribution to `log(log(y))`, evaluating the `1 - e` quantile with SciPy, exponentiating back to the original scale, and rounding to six decimals.

## Error Handling

Both public threshold functions raise `ValueError` for insufficient usable data or degenerate inputs. `bw_nrd0` raises `ValueError` when the finite data are constant.

## Dependencies

Required Python dependencies:

- `numpy`
- `scipy`

No runtime R dependencies remain.

## Validation

Validation is covered by:

- `tests/test_threshold_lib.py`
- golden output tests in `tests/test_regression_outputs.py`

For the current `mm`, `kp`, and `frogs` golden fixtures, non-Newick outputs match v1.3.9 exactly and Newick outputs are semantically identical.
