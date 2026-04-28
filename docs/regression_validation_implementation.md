# Regression Validation Implementation Note

This note records the formal validation added for TreeShrink v1.4.0.

## Implemented Tests

The test suite uses Python's standard-library `unittest` module and can be run with:

```bash
python -m unittest discover -s tests -p 'test_*.py'
```

Implemented files:

- `tests/test_threshold_lib.py`
- `tests/test_regression_outputs.py`
- `tests/test_runtime_logging.py`
- `tests/test_python3_compatibility.py`
- `tests/helpers/compare_newick_semantics.py`

## Threshold Validation

`tests/test_threshold_lib.py` validates the Python threshold implementation that replaced:

- `R_scripts/find_threshold_lkernel.R`
- `R_scripts/find_threshold_loglnorm.R`

The tests cover fixed numeric arrays, multiple threshold paths, insufficient data, and degenerate data.

## Golden Output Validation

`tests/test_regression_outputs.py` compares v1.4.0 output against tracked v1.3.9 golden fixtures for three datasets:

- `mm`
- `kp`
- `frogs`

Each dataset is tested in three modes:

- all-genes
- per-gene
- per-species

Fixture layout:

```text
tests/
  fixtures/
    inputs/
    golden/
  helpers/
    compare_newick_semantics.py
```

Comparison rules:

- removal-set `.txt` files are compared byte-for-byte
- `_summary.txt` files are compared byte-for-byte
- `.trees` files are compared as parsed Newick trees, not byte-for-byte
- generated runtime logs are ignored

The golden test was split into per-dataset/per-mode test methods so individual slow cases can be rerun directly.

## Runtime Logging Validation

`tests/test_runtime_logging.py` validates:

- `<prefix>.log` creation
- stdout/stderr tee behavior
- prefix collision handling
- `--force` overwrite behavior
- warning capture
- exception capture after logging starts
- lightweight `--version` and help behavior

## Python 3 Compatibility Validation

`tests/test_python3_compatibility.py` validates the TreeShrink-owned Python 3 fixes and confirms private vendored DendroPy import behavior.

## Compatibility Result

For the current `mm`, `kp`, and `frogs` fixtures:

- non-Newick outputs match v1.3.9 exactly
- Newick outputs are semantically identical
- runtime logging and Python 3 compatibility tests pass
