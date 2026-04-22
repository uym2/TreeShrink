# Future TreeSwift Migration

This note records future work to replace TreeShrink's remaining DendroPy usage with TreeSwift. This was not part of the v1.4.0 packaging cleanup.

## Current State

TreeShrink still uses vendored DendroPy 4.3.0 for tree parsing, tree writing, and some alignment-related paths. The vendored copy now lives privately under:

```text
treeshrink/_vendor/dendropy
```

TreeSwift is already used by `decompose.py` and `treeshrink.decompose_lib`.

## Migration Principles

- Preserve output behavior before changing implementation details.
- Replace one behavior surface at a time.
- Use golden output validation as the baseline.
- Compare Newick outputs semantically rather than byte-for-byte.
- Keep vendored DendroPy until all required behavior has been replaced and validated.

## Likely Migration Areas

- `run_treeshrink.py`
- `treeshrink/tree_lib.py`
- `treeshrink/Tree_extend.py`
- `treeshrink/filter_lib.py`
- `treeshrink/optimal_filter_lib.py`
- `treeshrink/alignment.py`
- related scripts that import TreeShrink tree/filter helpers

## Suggested Order

1. Inventory DendroPy usage by behavior.
2. Replace Newick parsing and writing first.
3. Replace traversal, pruning, and branch-length helpers.
4. Validate golden outputs after each slice.
5. Add focused tests for tree operations not covered by golden outputs.
6. Defer alignment-specific replacement until tree-only behavior is stable.
7. Remove vendored DendroPy after no runtime path imports it.

## Validation Baseline

Before removing DendroPy, run:

```bash
python -m unittest discover -s tests -p 'test_*.py'
```

Additional migration tests should cover branch length preservation, taxon labels with underscores, pruning behavior, tree writing, and rooted/unrooted assumptions.
