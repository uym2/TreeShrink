# Future TreeSwift Migration

This note records the longer-term plan to replace TreeShrink's DendroPy usage with TreeSwift in a future release. This migration is separate from the v1.4.0 Python 3 compatibility cleanup and the vendored DendroPy namespace move.

## Motivation

TreeShrink currently relies on a vendored DendroPy 4.3.0 copy for tree parsing, tree writing, and some alignment-related code paths. The vendored copy keeps current behavior stable, but it also carries compatibility and packaging costs.

Moving to TreeSwift would reduce dependency surface and align TreeShrink around a tree library already used by `decompose.py` and `treeshrink.decompose_lib`.

## Migration Principles

- Preserve TreeShrink output behavior before changing implementation details.
- Replace one behavior surface at a time.
- Use golden output validation as the compatibility baseline.
- Compare Newick outputs semantically rather than byte-for-byte during the transition.
- Keep the vendored DendroPy copy until all required behavior has been replaced and validated.

## Likely Migration Areas

TreeShrink-owned files that currently import or depend on DendroPy include:

- `run_treeshrink.py`
- `treeshrink/tree_lib.py`
- `treeshrink/Tree_extend.py`
- `treeshrink/filter_lib.py`
- `treeshrink/optimal_filter_lib.py`
- `treeshrink/alignment.py`
- related scripts that import TreeShrink tree/filter helpers

## Implementation Outline

1. Inventory DendroPy usage by behavior, not only by import.
2. Start with Newick parsing and writing, because this is central to TreeShrink outputs.
3. Replace tree traversal and pruning helpers with TreeSwift-backed equivalents.
4. Validate all existing golden output tests after each behavior slice.
5. Add focused tests for any tree operations not already covered by golden outputs.
6. Defer alignment-specific DendroPy replacement until tree-only behavior is stable.
7. Remove vendored DendroPy only after no TreeShrink runtime path imports it.

## Validation

The current regression suite provides the baseline:

- threshold unit tests
- runtime logging tests
- golden output tests for `mm`, `kp`, and `frogs`
- semantic Newick comparison for generated `.trees` files

Before removing DendroPy, run:

```bash
python -m unittest discover -s tests -p 'test_*.py'
```

For migration work, also add targeted tests around:

- branch length preservation
- taxon label preservation, including underscores
- pruning behavior
- tree writing behavior
- rooted and unrooted assumptions where relevant

## Release Boundary

This should not be part of the v1.4.0 packaging cleanup. The v1.4.0 path should keep vendored DendroPy under `treeshrink/_vendor/dendropy` and avoid external DendroPy dependencies. The TreeSwift migration should be a later release with its own validation cycle.
