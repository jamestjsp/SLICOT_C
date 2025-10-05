# agenet.md

This guide tailors the SLICOT repository workflow for the Droid coding agent.

## Interaction Principles
- Follow explicit user and system instructions exactly; avoid unsolicited suggestions.
- Keep replies concise (≤4 sentences when possible) and omit emojis unless requested.
- Use repository tools (`Read`, `ApplyPatch`, `TodoWrite`, etc.) as needed; maintain an up-to-date todo list for non-trivial tasks.
- Match existing code style, add only essential comments, and never expose sensitive information.

## Project Overview
SLICOT (Subroutine Library In COntrol Theory) provides:
- Fortran 77 routines in `src/`
- C wrappers in `src_c_wrapper/`
- GoogleTest suites in `tests/`
- Examples in `examples/`
- Modern CMake and legacy Makefile build systems

## Build System Commands
### CMake (preferred)
```bash
cmake --preset macos-x64-debug
cmake --preset macos-x64-release
cmake --preset linux-x64-debug
cmake --preset linux-x64-release
cmake --preset windows-x64-debug-intel
cmake --preset windows-x64-debug-mingw

cmake --build --preset macos-x64-debug-build
ctest --preset macos-x64-debug-test
cmake --build --preset macos-x64-debug-build --target install
```

### Legacy Makefiles
```bash
make -f makefile_Unix
nmake -f makefile
```

### Key CMake Options
- `SLICOT_BUILD_EXAMPLES`
- `SLICOT_BUILD_TESTS`
- `SLICOT_BUILD_C_WRAPPERS`
- `SLICOT_BUILD_SHARED_LIBS`

**Integer Size**: C wrappers rely on 32-bit `int`, supporting dimensions up to 2³¹ − 1.

**BLAS/LAPACK Threading**: Builds default to single-threaded MKL; keep computations single-threaded.

## Repository Architecture
- `src/`: Fortran routines (chapters A, B, D, F, I, M, N, S, T, U)
- `src_c_wrapper/`: C interoperability using `F77_FUNC`
- `src_aux/`: LAPACK compatibility helpers
- `include/`: Header files (`slicot_f77.h`, `slicot_utils.h`, etc.)
- `tests/`: GoogleTest suites with data fixtures

## C Wrapper Requirements
- Allocate/free workspace internally; do not expect caller-managed buffers.
- Support `row_major` flag by transposing through temporary column-major buffers when necessary.
- Derive Fortran leading dimensions from matrix row counts; compute C leading dimensions for row-major callers.
- Mirror Fortran validation: allow zero-dimension quick exits and tolerate `NULL` when Fortran ignores data.
- Use `SLICOT_EXPORT`, `CHECK_ALLOC`, and `goto cleanup` patterns consistently.

## Testing Guidelines
- Each wrapper requires a companion `*_test.cpp` using GoogleTest fixtures for row- and column-major paths.
- Embed small datasets directly; store larger datasets in `tests/data/` and load with `load_test_data_from_csv()`.
- Choose realistic tolerances (`1e-3` to `5e-3`) based on routine behavior and verify zero-dimension scenarios.
- Compare outputs to validated references or documentation where available.

## Workflow for Enhancements
1. Review Fortran docs/examples and `rep96-1.pdf` for parameter limits and workspace formulas.
2. Ensure validation matches Fortran expectations, especially for degenerate dimensions.
3. Place new Fortran code in `src/`, wrappers in `src_c_wrapper/`, headers in `include/`.
4. Update `cmake/slicot_sources.cmake` when adding sources.
5. Add unit tests and data as needed.
6. Run preset builds and relevant tests before completion unless user opts out.

## Utility References
- `slicot_transpose_to_fortran_with_ld()` / `slicot_transpose_to_c_with_ld()` for layout conversions.
- `load_test_data_from_csv()` for CSV fixtures.
- Macros: `CHECK_ALLOC`, `SLICOT_MEMORY_ERROR`, `MAX`, `MIN`.

## Additional Reminders
- Respect zero-sized array handling and packed/band storage conventions outlined in `rep96-1.pdf`.
- Investigate `tests/Known_Issues.txt` for acknowledged failures before reporting regressions.
- Before committing, inspect `git status` and `git diff --cached`; never push unless instructed.
- Do not modify documentation (README, etc.) unless explicitly requested by the user.
