# Packaging Implementation Note For TreeShrink v1.4.0

TreeShrink v1.4.0 removed the R/BMS runtime path, added Python threshold dependencies, and moved vendored DendroPy into TreeShrink's private namespace.

## Removed Legacy Files

- Removed `ez_setup.py`, the obsolete setuptools bootstrap script.
- Removed `dependencies/`, which only contained archived third-party tarballs:
  - `BMS_0.3.3.tar.gz`
  - `DendroPy-4.3.0.tar.gz`

`BMS_0.3.3.tar.gz` is no longer needed because v1.4.0 does not call R/BMS. `DendroPy-4.3.0.tar.gz` was redundant because the vendored source is included directly.

## Vendored DendroPy

Vendored DendroPy 4.3.0 was moved from the top-level `dendropy` package to:

```text
treeshrink/_vendor/dendropy
```

TreeShrink imports now use `treeshrink._vendor.dendropy`, and vendored DendroPy internal imports were updated to the private namespace.

Packaging checks confirmed:

- `find_packages()` includes `treeshrink._vendor.dendropy`
- `find_packages()` does not include top-level `dendropy`
- external `DendroPy` is not listed as a pip or conda dependency

## `setup.py`

`setup.py` now uses:

```python
packages=find_packages()
install_requires=[
    "treeswift",
    "numpy",
    "scipy",
]
python_requires=">=3.8"
```

`R_scripts/` remains in the source tree as historical/reference material, but it is no longer installed as a Python package or runtime package data.

## Conda Package

The conda recipe was updated for v1.4.0:

- version is `1.4.0`
- source revision is `v1.4.0`
- `r-base` was removed from host and run requirements
- `treeswift`, `numpy`, and `scipy` were added to run requirements
- BMS installation commands were removed from `build.sh` and `bld.bat`
- `run_treeshrink.py --version` was added as a conda command test

## Validation

Packaging checks completed:

```bash
python setup.py check
python -m pip install --no-deps --target /tmp/treeshrink_pkg_smoke .
```

The install smoke test confirmed that the installed package contains private vendored DendroPy and does not install top-level `dendropy`, `R_scripts`, or `dependencies`.

Full tests also passed:

```bash
python -m unittest discover -s tests -p 'test_*.py'
```

## Remaining

- Build the conda package during release preparation.
