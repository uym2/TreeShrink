# Runtime Logging Implementation Note

TreeShrink v1.4.0 automatically writes runtime terminal messages to a log file while still showing those messages on screen.

## Implemented Behavior

- Help-only and version-only invocations exit before logging starts, so they do not create output directories or log files.
- Normal runs write `<outdir>/<prefix>.log` beside the regular TreeShrink outputs.
- The log uses the final output prefix after collision handling.
- Standard output and standard error are mirrored into the same log file.
- Terminal output is preserved.

## Implementation Details

- `run_treeshrink.py` resolves the output directory and final prefix immediately after argument parsing and version handling.
- `prepare_output()` centralizes output directory creation and prefix collision handling.
- `TeeStream` mirrors writes to both the original terminal stream and the open log file.
- `start_runtime_logging()` replaces `sys.stdout` and `sys.stderr` before the launch messages are printed.
- `stop_runtime_logging()` restores both streams and closes the log file from the top-level `finally` block.

## Prefix Semantics

- If `--force` is not used and the requested prefix already exists, TreeShrink increments the prefix before opening the log.
- If `--force` is used, TreeShrink uses the requested prefix and overwrites the corresponding log.
- The log prefix matches the prefix used for summary, removal-set, tree, and alignment outputs.

## Validation Performed

- `python -m py_compile TreeShrink/run_treeshrink.py`
- Normal run on `test_data/mm10.trees` created `output.log`.
- Repeated run without `--force` created `output1.log`.
- Run with `--force` overwrote `output.log`.
- `--version` still printed only `1.4.0` without starting logging.
