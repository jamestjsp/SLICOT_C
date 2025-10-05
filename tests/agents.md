# tests/agents.md

## Purpose
- Implement GoogleTest suites for SLICOT C wrappers, covering both column- and row-major paths.
- Validate numerical outputs against verified references while honoring SLICOT quick-exit behaviors.

## Data Handling
- Prefer small, embedded datasets for fixtures; otherwise load CSV files from `tests/data/` using `load_test_data_from_csv()` and keep headers identical to the requested column lists.
- When parsing HTML "Program Data" / "Program Results" tables, read values sequentially as printed but verify layout against Fortran READ loops (`((..., I=1,rows), J=1,cols)`) to decide whether a transpose is required; reconcile discrepancies by consulting the companion `T*.f` driver when available.
- Apply `slicot_transpose_to_c_with_ld()` and related helpers to convert between storage orders, and retain the original column-major ordering when storing vectors for direct comparison with Fortran expectations.
- For datasets derived from `examples/*.dat`, benchmark archives, or external sources, cross-check dimension conventions and record layouts in `docs/rep96-1.pdf` Sections 2.2.6–2.2.10 and Chapter 4 before reshaping samples or creating CSV files.

## Test Coverage
- Assert wrapper validation for negative dimensions and confirm zero-dimension cases mirror underlying Fortran behavior (`INFO=0` quick exits where documented).
- Compare outputs with tolerances derived from routine sensitivity (start at `1e-3`, relax only as justified by repeated runs).
- Include scenario coverage for example problems, benchmark problems, and datasets referenced in `doc/*.html`, `examples/*.dat`, `examples/*.res`, and `benchmark_data/`, documenting any verified divergences from published results in test comments or commit notes rather than altering source data.

## Execution
- Build via `cmake --build --preset <platform>-debug-build` and execute `ctest --preset <platform>-debug-test` or run targeted binaries from `build/<preset>/tests/`.
- Keep diagnostics concise; prefer GTest assertions for validation instead of manual logging.
