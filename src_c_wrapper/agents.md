# src_c_wrapper/agents.md

## Purpose
- Implement C interfaces that invoke SLICOT Fortran routines through `F77_FUNC` while abstracting workspace management and layout conversion for callers.

## Core Requirements
- Allocate all workspace (`iwork`, `dwork`, etc.) internally; follow sizing rules from routine docs and `docs/rep96-1.pdf` Section 2.2.10. Support workspace queries when the Fortran routine permits.
- Only treat `DWORK(1)` (or similar slots) as the minimum required workspace when the routine documentation explicitly guarantees that behaviour; otherwise propagate `INFO=-23` without retrying.
- Honor the `row_major` flag: transpose inputs with `slicot_transpose_to_fortran_with_ld()` and restore outputs via `slicot_transpose_to_c_with_ld()` when `row_major == 1`.
- Compute Fortran leading dimensions from row counts regardless of caller layout while validating C-side LD arguments per `docs/CONTRIBUTING.md` guidelines.
- Align zero-dimension handling with the underlying routine (quick exits with `INFO=0` where documented) and accept `NULL` pointers only when the Fortran docs state the arrays are not referenced.

## Validation & Errors
- Perform argument checks before invoking Fortran; return negative indices for wrapper validation failures and propagate Fortran `INFO` otherwise.
- Use `CHECK_ALLOC` and `goto cleanup` patterns consistently; release every allocation before returning.
- Keep filenames lowercase, match existing style, and add comments only to clarify deviations from the common wrapper template.
