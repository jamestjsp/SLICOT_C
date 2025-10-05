# examples/agents.md

## Purpose
- Preserve official SLICOT example programs (`T*.f`) and datasets (`*.dat`, `*.res`) that demonstrate routine usage and expected outputs.

## Working With Example Problems
- Treat `.f` drivers and accompanying data as canonical references; do not rewrite unless maintainers approve.
- When extracting values for tests or wrappers, follow the READ order described in `docs/rep96-1.pdf` Chapter 4 and Appendix C to map records into matrices/vectors correctly.
- Document any verified deviations from published `.res` outputs in commit messages or dedicated test notes rather than editing the historical files.

## Data Guidance
- Before converting `.dat` files to CSV for `tests/data/`, confirm dimension ordering, sample stride, and complex-number conventions from the example’s HTML doc and `rep96-1` Sections 2.2.6–2.2.9.
- Maintain original formatting so other tools can replay the Fortran examples without modification.
