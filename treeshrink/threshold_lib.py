import numpy as np
from scipy.stats import norm


def bw_nrd0(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    n = x.size
    if n < 2:
        raise ValueError("Need at least 2 finite observations.")

    sd = np.std(x, ddof=1)
    q75, q25 = np.percentile(x, [75, 25], method="linear")
    iqr = q75 - q25
    lo = min(sd, iqr / 1.34)

    if lo <= 0:
        if sd > 0:
            lo = sd
        else:
            raise ValueError("Data are constant; bandwidth is zero.")

    return 0.9 * lo * n ** (-1 / 5)


def _r_bindist(x: np.ndarray, weights: np.ndarray, lo: float, up: float,
               n: int) -> np.ndarray:
    """Linear binning matching R's internal C_BinDist for in-range values."""
    y = np.zeros(2 * n, dtype=float)
    xdelta = (up - lo) / (n - 1)
    xpos = (x - lo) / xdelta
    ix = np.floor(xpos).astype(int)
    fx = xpos - ix

    in_left = (0 <= ix) & (ix < n)
    np.add.at(y, ix[in_left], (1.0 - fx[in_left]) * weights[in_left])

    ix_right = ix + 1
    in_right = (0 <= ix_right) & (ix_right < n)
    np.add.at(y, ix_right[in_right], fx[in_right] * weights[in_right])
    return y


def r_like_density_values(x: np.ndarray, adjust: float = 1.0, n: int = 512,
                          cut: float = 3.0) -> tuple[np.ndarray, np.ndarray, float]:
    """
    Match R stats::density.default for the default Gaussian/nrd0 path.

    The R implementation bins observations on an extended grid, convolves those
    bins with the Gaussian kernel via FFT, clamps tiny negative numerical noise,
    and interpolates back to the user-facing grid.
    """
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 2:
        raise ValueError("Need at least 2 finite observations.")

    n_user = n
    n = max(n, 512)
    if n > 512:
        n = 2 ** int(np.ceil(np.log2(n)))

    bw = bw_nrd0(x) * adjust
    from_ = np.min(x) - cut * bw
    to = np.max(x) + cut * bw
    lo = from_ - 4 * bw
    up = to + 4 * bw

    weights = np.full(x.size, 1.0 / x.size)
    y = _r_bindist(x, weights, lo, up, n)

    kords = np.linspace(0.0, 2.0 * (up - lo), 2 * n)
    kords[n + 1:] = -kords[n - 1:0:-1]
    kernel = norm.pdf(kords, scale=bw)

    conv = np.fft.ifft(np.fft.fft(y) * np.conj(np.fft.fft(kernel)))
    dens_ext = np.maximum(0.0, np.real(conv)[:n])

    xords = np.linspace(lo, up, n)
    grid = np.linspace(from_, to, n_user)
    dens = np.interp(grid, xords, dens_ext)

    return grid, dens, bw


def quantile_from_density_grid(grid: np.ndarray, dens: np.ndarray, p: float) -> float:
    """Match BMS::quantile.density for a single probability."""
    if not (0 <= p <= 1):
        raise ValueError("p must be in [0, 1].")

    grid = np.asarray(grid, dtype=float)
    dens = np.asarray(dens, dtype=float)
    if grid.ndim != 1 or dens.ndim != 1 or grid.size != dens.size:
        raise ValueError("grid and dens must be one-dimensional arrays of equal length.")
    if grid.size < 2:
        raise ValueError("Need at least 2 density grid points.")

    dx = grid[1] - grid[0]
    cdf = (np.cumsum(dens) - (dens - dens[0]) / 2.0) * dx
    total = cdf[-1]
    if total <= 0 or not np.isfinite(total):
        raise ValueError("Density area must be positive and finite.")
    cdf = cdf / total

    iii = int(np.sum(cdf <= p))
    if iii == cdf.size:
        return float("inf")
    if iii == 0:
        return float("-inf")

    left = iii - 1
    right = iii
    return float(
        grid[right]
        + ((cdf[right] - p) / (cdf[right] - cdf[left]))
        * (grid[left] - grid[right])
    )


def threshold_l_kernel(y, e=0.05):
    y = np.asarray(y, dtype=float)
    x = y[(y > 0) & np.isfinite(y)]
    if x.size < 2:
        raise ValueError("Need at least 2 positive finite values.")

    logx = np.log(x)
    grid, dens, _ = r_like_density_values(logx, adjust=1.0, n=512, cut=3.0)
    q = quantile_from_density_grid(grid, dens, 1 - e)
    return round(float(np.exp(q)), 6)


def threshold_loglnorm(y, e=0.05):
    """
    Match R_scripts/find_threshold_loglnorm.R.

    The R script computes:
        exp(qlnorm(1 - e,
                   meanlog = mean(log(log(y[y > 1]))),
                   sdlog = sd(log(log(y[y > 1])))))
    """
    y = np.asarray(y, dtype=float)
    logy = np.log(y)
    x = logy[(logy > 0) & np.isfinite(logy)]
    if x.size < 2:
        raise ValueError("Need at least 2 values with positive finite logs.")

    logx = np.log(x)
    meanlog = np.mean(logx)
    sdlog = np.std(logx, ddof=1)
    q = np.exp(meanlog + sdlog * norm.ppf(1 - e))
    return round(float(np.exp(q)), 6)
