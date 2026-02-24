# SUPPA2 splicing workflow (events → filtering → PSI → differential splicing)

This repository documents a reproducible workflow to run SUPPA2 on RNA-seq transcript quantification (TPM) data:
1) generate alternative splicing events from a GTF,
2) validate and merge TPM files across replicates,
3) filter events by transcript availability across all replicates,
4) compute PSI per event, and
5) perform differential splicing analyses across timepoints.

## Documentation
- Full step-by-step pipeline: **[`docs/pipeline.md`](docs/pipeline.md)**

## Requirements
- Python >= 3.10
- SUPPA2 (tested with v2.4)
- Python packages: pandas

## Repository structure
- `scripts/` — helper scripts (QC, filtering, merging)
- `docs/` — detailed documentation
- `examples/` — optional toy inputs/outputs

## Quickstart (outline)
```bash
# 1) Create and activate a virtual environment
python3 -m venv .venv
source .venv/bin/activate

# 2) Install Python dependencies
pip install -r requirements.txt

# 3) Run the pipeline following docs/pipeline.md
```
