# SUPPA2 splicing workflow (events → filtering → PSI → differential splicing)

This repository provides a reproducible, documented workflow to run **SUPPA2** on transcript quantification (TPM) data:
- generate alternative splicing events from a GTF,
- validate and merge TPM files across replicates,
- filter events by transcript availability across all replicates,
- compute PSI per event,
- perform differential splicing analyses across timepoints.

## Documentation
- Full step-by-step pipeline: **[`docs/pipeline.md`](docs/pipeline.md)**

## Requirements
- Python >= 3.10
- SUPPA2 (tested with v2.4)
- Python packages: pandas

## Repository structure
- `scripts/` — helper scripts (QC, filtering, merging)
- `docs/` — detailed documentation

## Quickstart (outline)
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# Follow the full guide:
# docs/pipeline.md
```
## Citation
If you use this workflow, please cite SUPPA2:
https://github.com/comprna/SUPPA

## License
MIT License. See [`LICENSE`](LICENSE).
