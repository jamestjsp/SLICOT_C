# benchmark_data/agents.md

## Purpose
- Store datasets used for performance and regression benchmarks of SLICOT routines.

## Data Stewardship
- Keep raw benchmark inputs untouched; add new datasets in clearly labeled subdirectories with metadata identifying source problems.
- When interpreting files, consult `docs/rep96-1.pdf` Section 2.2 (dimension conventions) and Chapter 4 for record layouts to ensure consistency with example and test ingestion pipelines.
- Note whether data originate from example problems, industrial cases, or synthetic generators so downstream tests can cite provenance.

## Usage Guidance
- Coordinate with `tests/` maintainers before converting benchmark data into CSV or alternative formats to avoid duplication and ensure compatible header naming.
- Document any preprocessing scripts or scaling applied to benchmark matrices in accompanying notes or commit messages.
