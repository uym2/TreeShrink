# Regression Validation Implementation Note

This note records the validation approach used while moving TreeShrink v1.4.0 away from R-based thresholding and adding automatic runtime logs. It also defines the regression checks future formal tests should preserve.

## Validation Goals

- Protect threshold behavior after replacing R scripts with `treeshrink.threshold_lib`.
- Confirm TreeShrink output compatibility across future versions.
- Verify runtime log creation without changing terminal behavior.
- Avoid brittle byte-level comparisons for Newick strings that represent the same tree.

## Threshold Validation

The Python threshold implementation replaces two R scripts:

- `R_scripts/find_threshold_lkernel.R`
- `R_scripts/find_threshold_loglnorm.R`

Formal unit tests should cover:

- `threshold_l_kernel(y, e=0.05)` on fixed numeric arrays
- `threshold_loglnorm(y, e=0.05)` on fixed numeric arrays
- multiple quantile values where useful
- insufficient usable data raising `ValueError`
- constant or degenerate data raising `ValueError`

The threshold tests currently use Python's standard-library `unittest` module so they can run without adding a new test dependency. Future integration tests can continue with `unittest` or move to `pytest` if the project adopts it.

Current test command:

```bash
python -m unittest discover -s tests -p 'test_*.py'
```

Current implemented suite:

- `tests/test_threshold_lib.py`
- `tests/test_regression_outputs.py`
- `tests/test_runtime_logging.py`

## Golden Output Validation

The golden output suite covers three tracked input datasets, `mm`, `kp`, and `frogs`, across three TreeShrink modes:

- all-genes
- per-gene
- per-species

The validation compares generated output directories for v1.3.9 and v1.4.0.

Golden tests are self-contained and use tracked fixtures. Input trees and v1.3.9 golden outputs are stored under `tests/fixtures/`.

Implemented test layout:

```text
tests/
  fixtures/
    inputs/
      mm.trees
      kp.trees
      frogs.trees
    golden/
      mm_allgenes_v139/
      mm_pergene_v139/
      mm_perspecies_v139/
      kp_allgenes_v139/
      kp_pergene_v139/
      kp_perspecies_v139/
      frogs_allgenes_v139/
      frogs_pergene_v139/
      frogs_perspecies_v139/
  helpers/
    compare_newick_semantics.py
  test_regression_outputs.py
```

For each dataset and mode, `tests/test_regression_outputs.py` runs TreeShrink against the matching tracked input tree in a temporary output directory and compares against the matching v1.3.9 golden directory:

- removal-set `.txt` files byte-for-byte
- `_summary.txt` files byte-for-byte
- `.trees` files semantically as parsed Newick trees, not byte-for-byte

Generated runtime logs should be ignored by golden output comparison because v1.3.9 did not produce log files.

Acceptance criteria:

- all non-Newick outputs match exactly
- all Newick outputs represent semantically identical trees with matching branch lengths within a small tolerance

Implementation notes:

- run generated outputs under temporary directories
- keep golden expected files small enough for routine CI
- keep semantic Newick comparison logic in `tests/helpers/`
- compare command exit codes and captured stdout/stderr in addition to output files

Compatibility baseline:

- for the current `mm`, `kp`, and `frogs` outputs, non-Newick outputs match exactly
- all parsed Newick trees are semantically identical
- future golden tests should treat this behavior as the v1.4.0 compatibility baseline

## Runtime Logging Validation

Runtime logging is covered by `tests/test_runtime_logging.py`, using a small self-contained Newick fixture under temporary output directories.

The implemented tests assert:

- `<prefix>.log` is created
- the log contains launch and called-as messages
- warnings, such as low tree-count mode switching, appear in the log
- the final output-location message appears in the log
- terminal output is still produced

Prefix behavior to preserve:

- first run creates `output.log`
- second run without `--force` creates `output1.log`
- run with `--force` overwrites `output.log`

Special CLI behavior should remain lightweight:

- `--version` prints only the version and creates no log
- help-only invocation creates no output directory and no log
- exceptions after logging starts are written to both stderr and `<prefix>.log`
