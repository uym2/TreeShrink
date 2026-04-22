# Python 3 Compatibility Implementation Note

This note covers the immediate Python 3 compatibility cleanup for TreeShrink-owned code in v1.4.0. It intentionally does not refactor or replace the vendored `dendropy` package in this step.

## Scope

The current implementation should fix modern-Python hazards in TreeShrink's own runtime code while preserving behavior and keeping changes small.

Included:

- replace invalid universal-newline file mode in `treeshrink/alignment.py`
- replace Python 2 file iteration in `treeshrink/sequence_lib.py`
- add focused tests for the affected code paths

Not included in the first small compatibility patch:

- broad modernization of vendored `dendropy`
- replacing DendroPy usage with TreeSwift

## DendroPy Vendor Policy

For v1.4.0, TreeShrink continues using the vendored `dendropy` package included in the repository, but the vendored code has moved out of the top-level package namespace.

This keeps TreeShrink insulated from upstream DendroPy API and behavior changes. It also means packaging should not declare external `DendroPy` as an installation dependency while the vendored copy remains active.

Policy and implementation:

- keep the vendored DendroPy code
- move it from top-level `dendropy` to `treeshrink/_vendor/dendropy`
- do not add `DendroPy` to `setup.py` `install_requires`
- do not add external `dendropy` to the conda runtime requirements
- update TreeShrink imports to use `treeshrink._vendor.dendropy`
- update vendored DendroPy internal imports to use the private namespace

## Compatibility Scan Result

A scan of TreeShrink-owned Python files found two actionable runtime issues.

### `treeshrink/alignment.py`

Problem:

```python
fileobj = open(seq_fn, 'rU')
```

The `'rU'` universal-newline mode is invalid in modern Python. Python 3 already handles universal newlines in normal text mode, so the behavior-preserving replacement is:

```python
with open(seq_fn, 'r') as fileobj:
    sd.read(
        fileobj,
        file_format=file_format,
        datatype=datatype,
        filename=seq_fn,
        careful_parse=careful_parse,
    )
```

Implementation notes:

- keep `file_format`, `datatype`, `filename`, and `careful_parse` handling unchanged
- use a context manager so the file is closed even when parsing fails
- do not change missing-data validation or exception behavior
- avoid larger refactors in `MultiLocusDataset.read_files`

### `treeshrink/sequence_lib.py`

Problem:

```python
taxon_dict[line[1:-1]] = gap_rm(f.next().rstrip())
```

Python 3 file objects do not expose `.next()`. The direct replacement is:

```python
taxon_dict[line[1:-1]] = gap_rm(next(f).rstrip())
```

Implementation notes:

- preserve current FASTA parsing assumptions
- preserve existing gap removal behavior
- avoid changing duplicate-taxon behavior or multiline FASTA handling in this cleanup
- use a context manager only if it can be done without changing parsing semantics

## Tests

Add focused tests under `tests/` for the two compatibility fixes.

Recommended test cases:

- `hash_taxon_seq` reads a simple FASTA file and removes gaps from sequences
- `MultiLocusDataset.read_files` can read a FASTA file under modern Python without using `'rU'`

The tests should use temporary files and standard-library `unittest`, matching the existing test suite.

Suggested verification command:

```bash
python -m unittest discover -s tests -p 'test_*.py'
```

## Vendored DendroPy Limitation

The repository includes a vendored `dendropy` package. The vendored copy reports version `4.3.0`, which is much older than current DendroPy releases.

The compatibility scan found legacy constructs inside vendored `dendropy`, including:

- `xrange`
- `open(..., "rU")`
- `collections.Mapping`
- Python 2 compatibility helpers such as `unicode`

Basic TreeShrink usage currently works with this vendored copy, and the regression tests exercise the TreeShrink paths that use DendroPy for Newick parsing. However, less-used DendroPy modules may behave incorrectly under modern Python.

Implementation note:

- patch only vendored DendroPy paths TreeShrink actually uses if urgent modern-Python breakages appear during the namespace migration
- otherwise avoid broad behavioral edits to vendored DendroPy
- remove the vendored copy only after replacing TreeShrink's DendroPy usage in a later release

## Vendored DendroPy Namespace Migration

The previous top-level vendored `dendropy` package could conflict with a user's installed DendroPy. It has been moved into TreeShrink's private namespace.

Target layout:

```text
treeshrink/
  _vendor/
    __init__.py
    dendropy/
      __init__.py
      ...
```

Implementation completed:

- moved the current top-level `dendropy/` directory to `treeshrink/_vendor/dendropy/`
- added `treeshrink/_vendor/__init__.py`
- updated TreeShrink-owned imports from `dendropy...` to `treeshrink._vendor.dendropy...`
- updated vendored DendroPy internal absolute imports from `dendropy...` to `treeshrink._vendor.dendropy...`
- added a test confirming TreeShrink imports the private vendored package
- confirmed `find_packages()` no longer includes a top-level `dendropy` package
- reran the golden regression suite after the import migration

Remaining packaging smoke test:

- build/install the package in a clean environment and confirm installing TreeShrink does not affect any separately installed external `dendropy`

This migration touched many imports and has a larger packaging surface than the small Python 3 compatibility fixes.

## Related Future Work

Replacing DendroPy usage with TreeSwift is intentionally separate from this compatibility cleanup. See `docs/treeswift_migration.md` for the migration note.
