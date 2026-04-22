# Packaging Notes For TreeShrink v1.4.0

This note records packaging cleanup identified after the v1.4.0 runtime and validation changes. The main packaging impact is that v1.4.0 replaces R-based threshold estimation with `treeshrink.threshold_lib`, which introduces Python dependencies and removes the runtime need for R and the R `BMS` package.

## `ez_setup.py`

`ez_setup.py` is an old vendored setuptools bootstrap script. It checks for setuptools, downloads ancient setuptools egg files when needed, and can be run directly to install or upgrade setuptools/EasyInstall.

The file appears redundant for current TreeShrink packaging:

- no tracked file imports `ez_setup`
- `setup.py` imports setuptools directly with `from setuptools import setup, find_packages`
- README installation instructions call `python setup.py install`
- conda packaging invokes `setup.py` directly
- no `setup.cfg`, `pyproject.toml`, or `MANIFEST.in` gives `ez_setup.py` a special role

Decision:

- remove `ez_setup.py`

The project no longer needs to preserve the old setuptools bootstrap workflow.

## `dependencies/`

The `dependencies/` directory currently contains archived third-party packages:

- `BMS_0.3.3.tar.gz`
- `DendroPy-4.3.0.tar.gz`

Decision:

- remove the `dependencies/` directory

Rationale:

- `BMS_0.3.3.tar.gz` is obsolete because v1.4.0 no longer uses R/BMS at runtime
- `DendroPy-4.3.0.tar.gz` is redundant because the vendored DendroPy source is already present in the repository and will move to `treeshrink/_vendor/dendropy`
- `_vendor` should contain importable vendored runtime code, not old source tarballs
- conda build scripts should stop referencing `dependencies/BMS_0.3.3.tar.gz`

## Vendored DendroPy Migration

TreeShrink currently includes a top-level vendored `dendropy` package. The vendored copy reports version `4.3.0`. This was intentional: TreeShrink depends on DendroPy behavior that may shift across upstream releases.

For v1.4.0, TreeShrink keeps using this vendored DendroPy copy, but it has moved from the top-level package name `dendropy` to TreeShrink's private namespace:

```text
treeshrink/_vendor/dendropy
```

Do not add external `DendroPy` to `install_requires`.

Packaging result:

- `find_packages()` includes `treeshrink._vendor.dendropy`
- `find_packages()` no longer includes a top-level `dendropy`
- TreeShrink should install its pinned DendroPy copy as a private package
- `setup.py` should not ask pip to install a second, newer DendroPy
- conda runtime requirements should not include external `dendropy`
- TreeShrink should no longer install a top-level package named `dendropy`

Reason:

- keeping the old top-level vendored package can shadow or conflict with a user's separately installed DendroPy
- moving it under `treeshrink._vendor` keeps TreeShrink pinned to known DendroPy behavior without polluting the global package namespace

Migration completed:

- moved the vendored code under `treeshrink/_vendor/dendropy`
- added `treeshrink/_vendor/__init__.py`
- updated TreeShrink imports to use the vendored namespace
- updated vendored DendroPy internal absolute imports from `dendropy...` to `treeshrink._vendor.dendropy...`
- stopped installing a top-level `dendropy` package from TreeShrink
- added a test proving TreeShrink imports the private vendored DendroPy copy

## `setup.py`

The current `setup.py` does not declare runtime dependencies. For v1.4.0, it should declare the external Python packages used by the installed scripts and `treeshrink` modules.

Recommended runtime dependencies:

```python
install_requires=[
    "treeswift",
    "numpy",
    "scipy",
]
```

Dependency rationale:

- DendroPy is intentionally omitted because TreeShrink vendors DendroPy 4.3.0 privately under `treeshrink._vendor`
- `treeswift`: `decompose.py` and `treeshrink.decompose_lib`
- `numpy`: threshold calculations and runtime data loading
- `scipy`: `scipy.stats.norm` in `treeshrink.threshold_lib`

Recommended metadata update:

```python
python_requires=">=3.8",
```

Use `>=3.6` instead only if old Python support is still a project goal.

The current package list still includes `R_scripts`:

```python
packages=find_packages() + ["R_scripts"]
package_data={"": recursive_list_dir("R_scripts")}
```

For v1.4.0, these R scripts are no longer used at runtime. Recommended cleanup:

```python
packages=find_packages()
```

Keep `R_scripts` packaged only if the project wants to distribute them as historical or reference files.

## Conda Package

The conda recipe under `conda_package/treeshrink/` still reflects v1.3.9 packaging.

Current stale items:

- `version: v1.3.9`
- `git_rev: master`
- `r-base >=4.0` in both host and run requirements
- `build.sh` and `bld.bat` install `dependencies/BMS_0.3.3.tar.gz` with `R CMD INSTALL`

Recommended `meta.yaml` shape:

```yaml
package:
  name: treeshrink
  version: 1.4.0

source:
  git_url: https://github.com/uym2/TreeShrink.git
  git_rev: v1.4.0

requirements:
  host:
    - python >=3.8
    - setuptools
  run:
    - python >=3.8
    - treeswift
    - numpy
    - scipy

test:
  imports:
    - treeshrink
  commands:
    - run_treeshrink.py --version

about:
  home: https://github.com/uym2/TreeShrink
  license: GPL-3.0-or-later
  license_file: LICENSE
  summary: Fast and accurate detection of outlier long branches in phylogenetic trees
```

If building locally during release preparation, use this source block instead:

```yaml
source:
  path: ../../
```

Recommended `build.sh`:

```bash
$PYTHON setup.py install --single-version-externally-managed --record=record.txt
```

Recommended `bld.bat`:

```bat
"%PYTHON%" setup.py install --single-version-externally-managed --record=record.txt
if errorlevel 1 exit 1
```

## Follow-Up Checklist

- remove `ez_setup.py`
- remove `dependencies/`
- add Python dependencies to `setup.py`
- decide whether `R_scripts` should still be included in source/package artifacts
- update conda recipe version and source revision for v1.4.0
- remove R/BMS install steps from conda build scripts
- add conda runtime dependencies: `treeswift`, `numpy`, `scipy`
- keep external `dendropy` out of pip and conda dependencies while the vendored copy is used
- move vendored DendroPy from top-level `dendropy` to `treeshrink/_vendor/dendropy`
- run package install smoke tests after updating packaging files
