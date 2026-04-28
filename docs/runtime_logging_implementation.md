# Runtime Logging Implementation Note

TreeShrink v1.4.0 writes runtime terminal messages to a log file while still showing them on screen.

## Implemented Behavior

- Help-only and version-only invocations exit before logging starts.
- Normal runs write `<outdir>/<prefix>.log`.
- The log uses the final output prefix after collision handling.
- Standard output and standard error are mirrored into the same log file.
- Terminal output is preserved.

## Implementation

`run_treeshrink.py` now:

- resolves output directory and final prefix immediately after argument parsing
- centralizes output directory creation and prefix collision handling in `prepare_output()`
- mirrors streams with `TeeStream`
- starts logging before launch messages are printed
- restores streams and closes the log in a top-level `finally` block
- prints tracebacks to both stderr and the log for exceptions after logging starts

## Prefix Semantics

- Without `--force`, existing output prefixes are auto-incremented before opening the log.
- With `--force`, the requested prefix is used and the log is overwritten.
- The log prefix matches summary, removal-set, tree, and alignment outputs.

## Validation

`tests/test_runtime_logging.py` covers log creation, prefix collision behavior, `--force`, warning capture, exception capture, and early `--version`/help exits.
