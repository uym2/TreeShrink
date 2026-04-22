# Threshold Library Implementation Note

`treeshrink/threshold_lib.py` replaces the R-based threshold calculations used by earlier TreeShrink versions. The goal is to keep v1.4.0 behavior compatible with the previous R scripts while removing the runtime dependency on `Rscript` and the R `BMS` package.

## Implemented Integration

`run_treeshrink.py` now imports:

- `threshold_l_kernel` for per-species and all-genes thresholding.
- `threshold_loglnorm` for per-gene thresholding.

These replace the previous subprocess calls to:

- `R_scripts/find_threshold_lkernel.R`
- `R_scripts/find_threshold_loglnorm.R`

## Kernel Threshold Implementation

`threshold_l_kernel` implements the previous R expression:

```r
x = y[y > 0]
exp(quantile.density(density(log(x), adjust = 1), p = 1 - e))
```

The Python implementation mirrors the default R `stats::density` Gaussian kernel path closely:

- filters to positive finite input values
- applies `log`
- computes R's `bw.nrd0` bandwidth rule
- performs R-like linear binning on the extended density grid
- applies Gaussian convolution with FFT
- interpolates back to the user-facing density grid
- computes the density quantile using the same trapezoid-style cumulative-density interpolation as `BMS::quantile.density`
- exponentiates and rounds the final threshold to six decimals

This preserves the old kernel-density threshold behavior without loading R or `BMS`.

## Log-Lognormal Threshold Implementation

`threshold_loglnorm` implements the previous R script:

```r
x = y[y > 0]
threshold = qlnorm(p = 1 - e, sdlog = sd(log(x)), meanlog = mean(log(x)))
exp(threshold(log(d$V1), e = e))
```

After simplifying the nested call, the Python implementation:

- takes `log(y)`
- keeps values whose log is positive and finite
- fits a lognormal distribution to `log(log(y))`
- evaluates the `1 - e` quantile using SciPy's normal quantile
- exponentiates back to the original scale
- rounds the final threshold to six decimals

## Implemented Error Handling

Both public threshold functions validate that there are at least two usable finite observations. They raise `ValueError` for insufficient or degenerate input rather than silently returning an invalid threshold.

`bw_nrd0` also raises `ValueError` when the finite data are constant, because R's density bandwidth would collapse to zero.

## Runtime Dependencies

The implementation depends on:

- `numpy`
- `scipy.stats.norm`

It intentionally does not depend on:

- `Rscript`
- R `stats`
- R `BMS`

## Validation Performed

Compatibility with older TreeShrink behavior was checked by comparing v1.3.9 and v1.4.0 outputs:

- compare v1.3.9 and v1.4.0 `.txt` and summary outputs byte-for-byte
- compare `.trees` outputs semantically, because equivalent Newick trees can be serialized differently
- use the repository comparison helpers:

```bash
/home/uym2/compare_treeshrink_outputs.sh
```

For the current Mammals test outputs, the non-Newick outputs match exactly and the Newick trees are semantically identical.
