# Python 3 Compatibility Implementation Note

This note records the Python 3 cleanup completed for TreeShrink-owned code in v1.4.0.

## TreeShrink-Owned Fixes

Two Python 2-era runtime hazards were fixed.

In `treeshrink/alignment.py`, invalid universal-newline mode:

```python
open(seq_fn, "rU")
```

was replaced with normal text mode and a context manager:

```python
with open(seq_fn, "r") as fileobj:
    ...
```

In `treeshrink/sequence_lib.py`, Python 2 file iteration:

```python
f.next()
```

was replaced with:

```python
next(f)
```

These changes preserved existing parsing behavior while making the affected paths work on modern Python.

## Vendored DendroPy Namespace

Vendored DendroPy remains in the repository for behavioral stability, but it no longer installs as a top-level `dendropy` package. It now lives under:

```text
treeshrink/_vendor/dendropy
```

TreeShrink and vendored DendroPy imports were updated to use `treeshrink._vendor.dendropy`.

## Tests

`tests/test_python3_compatibility.py` covers:

- `hash_taxon_seq()` reading FASTA under Python 3
- `MultiLocusDataset.read_files()` reading FASTA without `rU`
- TreeShrink importing private vendored DendroPy

The full test suite passed:

```bash
python -m unittest discover -s tests -p 'test_*.py'
```

## Remaining

Vendored DendroPy 4.3.0 still contains older compatibility patterns in less-used modules, such as `xrange`, `open(..., "rU")`, and `collections.Mapping`. The exercised TreeShrink paths pass. Broader DendroPy modernization is deferred because DendroPy will eventually be replaced by TreeSwift; see `docs/treeswift_migration.md`.
