# src/agents.md

## Purpose
- Maintain original SLICOT Fortran 77 routines organized by chapter (A–U) without altering established naming or calling conventions.
- Ensure any supplemental Fortran additions follow `rep96-1.pdf` interface, dimension, and workspace requirements.

## Implementation Notes
- Review the routine’s HTML documentation in `doc/` and cross-reference `docs/rep96-1.pdf` Sections 2.2–2.3 for argument order, leading-dimension rules, and zero-dimension handling.
- Preserve external entry points exactly; new code must keep arguments in column-major order and respect declared workspaces.
- Comment sparingly—match existing SLICOT style and rely on inline documentation only when clarifying deviations authorized by maintainers.

## Data & Examples
- When routines depend on example or benchmark datasets, confirm the Fortran READ patterns against Chapter 4 and Appendix C of `rep96-1.pdf` before modifying data flows.
- Align test expectations with corresponding `.res` files and flag discrepancies for downstream wrapper/test authors instead of altering reference outputs silently.
