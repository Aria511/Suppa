# SUPPA2 pipeline documentation

This document describes a step-by-step workflow for alternative splicing analysis with **SUPPA2** using:
- a gene annotation file (GTF),
- transcript expression quantification files (TPM),
- and replicate/timepoint comparisons.

> Note: paths in this document use placeholders. Replace them with your local paths.

---

## 0) Environment setup

### Check Python version
```bash
python --version
```

### Create and activate a virtual environment
```bash
python3 -m venv .venv
source .venv/bin/activate
```

### Install Python dependencies
```bash
pip install -r requirements.txt
```

---

## 1) SUPPA2 installation / usage

Download SUPPA2 v2.4:
https://github.com/comprna/SUPPA/releases/tag/v2.4

If you want to run SUPPA2 as `suppa` from anywhere, you can either:
- call it explicitly: `python /path/to/SUPPA-2.4/suppa.py ...`, or
- define a shell alias (example):

```bash
alias suppa='python /path/to/SUPPA-2.4/suppa.py'
```

Check that it works:
```bash
suppa -h
```

---

## 2) Define project paths (recommended)

Set a few environment variables to avoid hard-coded paths:

```bash
PROJECT_DIR=/path/to/project
GTF=$PROJECT_DIR/annotation.gtf

TPM_DIR=$PROJECT_DIR/tpm_files
EVENTS_DIR=$PROJECT_DIR/events
FILTERED_DIR=$PROJECT_DIR/filtered
PSI_DIR=$PROJECT_DIR/psi
DIFF_DIR=$PROJECT_DIR/diffsplice

mkdir -p "$EVENTS_DIR" "$FILTERED_DIR" "$PSI_DIR" "$DIFF_DIR"
```

---

## Inputs

Minimum required inputs:
- **GTF annotation**: `annotation.gtf`
- **TPM quantifications**: sample-level `.tpm` files (two-column format) or joined matrices produced by SUPPA2 `joinFiles`

Expected TPM file naming (example):
- Replicate group A: `A1.tpm`, `A2.tpm`, `A3.tpm`, `A4.tpm`
- Replicate group B: `B1.tpm`, `B2.tpm`, `B3.tpm`, `B4.tpm`
- Replicate group C: `C1.tpm`, `C2.tpm`, `C3.tpm`, `C4.tpm`

---

## Outputs (main)

Key outputs produced by the workflow:
- Event files: `events.ioe_*` and optionally a merged `events.ioe`
- Filtered events: `filtered_events_all_replicates.ioe`
- PSI matrices: `A_psi.psi`, `B_psi.psi`, `C_psi.psi`
- Timepoint-merged PSI/TPM: `time_*.psi`, `time_*.tpm`
- Differential splicing results: `diffsplice/*` (e.g., `*_dPSI`, `*_psivec`)

## 3) Generate alternative splicing events

Generate events from the GTF:

```bash
suppa generateEvents \
  -i "$GTF" \
  -o "$EVENTS_DIR/events" \
  -f ioe \
  -e SE SS MX RI FL
```

SUPPA2 outputs one `.ioe` file per event type.

### (Optional) Merge multiple `.ioe` files into a single `events.ioe`
If you need a single combined file:

```bash
first_file=$(ls "$EVENTS_DIR"/events.ioe_*.ioe | head -n 1)
head -n 1 "$first_file" > "$EVENTS_DIR/events.ioe"

for file in "$EVENTS_DIR"/events.ioe_*.ioe; do
  tail -n +2 "$file" >> "$EVENTS_DIR/events.ioe"
done

```

The resulting file has columns like:
- `seqname`
- `gene_id`
- `event_id`
- `alternative_transcripts`
- `total_transcripts`

---

## 4) TPM files: validation and merging replicates

### 4.1 Validate TPM format and numeric values
Before merging TPM files, ensure they:
- contain numeric expression values
- use consistent columns

Run the helper script (see `scripts/check_tpm_values.py`):

```bash
python scripts/check_tpm_values.py --input-dir "$TPM_DIR" --pattern "*.tpm"
```

### 4.2 Merge replicate TPM files
Example (4 samples per replicate group):

```bash
# Replicate group A
suppa joinFiles -f tpm -i "$TPM_DIR/A1.tpm" "$TPM_DIR/A2.tpm" "$TPM_DIR/A3.tpm" "$TPM_DIR/A4.tpm" -o "$TPM_DIR/A_all"

# Replicate group B
suppa joinFiles -f tpm -i "$TPM_DIR/B1.tpm" "$TPM_DIR/B2.tpm" "$TPM_DIR/B3.tpm" "$TPM_DIR/B4.tpm" -o "$TPM_DIR/B_all"

# Replicate group C
suppa joinFiles -f tpm -i "$TPM_DIR/C1.tpm" "$TPM_DIR/C2.tpm" "$TPM_DIR/C3.tpm" "$TPM_DIR/C4.tpm" -o "$TPM_DIR/C_all"
```

Merged TPM structure (example):
```
sample1 sample2 sample3 sample4
transcript1 <expr> <expr> <expr> <expr>
...
```

---

## 5) Filter events by transcripts shared across all replicates

Goal: keep only events whose transcripts are present in **all** replicate expression matrices.

Run the helper script (see `scripts/analyze_all_replicates.py`):

```bash
python scripts/analyze_all_replicates.py \
  --events-ioe "$EVENTS_DIR/events.ioe" \
  --tpm-a "$TPM_DIR/A_all.tpm" \
  --tpm-b "$TPM_DIR/B_all.tpm" \
  --tpm-c "$TPM_DIR/C_all.tpm" \
  --out-dir "$FILTERED_DIR"
```

Output:
- `filtered_events_all_replicates.ioe`

---

## 6) PSI quantification per event

Compute PSI per event for each replicate group:

```bash
# A
suppa psiPerEvent \
  -i "$FILTERED_DIR/filtered_events_all_replicates.ioe" \
  -e "$TPM_DIR/A_all.tpm" \
  -o "$PSI_DIR/A_psi"

# B
suppa psiPerEvent \
  -i "$FILTERED_DIR/filtered_events_all_replicates.ioe" \
  -e "$TPM_DIR/B_all.tpm" \
  -o "$PSI_DIR/B_psi"

# C
suppa psiPerEvent \
  -i "$FILTERED_DIR/filtered_events_all_replicates.ioe" \
  -e "$TPM_DIR/C_all.tpm" \
  -o "$PSI_DIR/C_psi"
```

PSI values are in `[0, 1]`:
- `0` = alternative isoform not used
- `1` = alternative isoform fully used

---

## 7) Merge PSI/TPM by timepoint (for differential analyses)

If your design has multiple timepoints (e.g., 0h, 16h, 20h, 24h) and A/B/C represent biological replicates, merge them per timepoint.

Run the helper script (see `scripts/merge_replicates.py`):

```bash
python scripts/merge_replicates.py \
  --psi-a "$PSI_DIR/A_psi.psi" \
  --psi-b "$PSI_DIR/B_psi.psi" \
  --psi-c "$PSI_DIR/C_psi.psi" \
  --tpm-a "$TPM_DIR/A_all.tpm" \
  --tpm-b "$TPM_DIR/B_all.tpm" \
  --tpm-c "$TPM_DIR/C_all.tpm" \
  --out-psi-dir "$PSI_DIR" \
  --out-tpm-dir "$TPM_DIR"
```

Expected outputs:
- PSI: `time_0h.psi`, `time_16h.psi`, `time_20h.psi`, `time_24h.psi`
- TPM: `time_0h.tpm`, `time_16h.tpm`, `time_20h.tpm`, `time_24h.tpm`

---

## 8) Differential splicing analysis (timepoint comparisons)

Example comparisons:

```bash
# 0h vs 16h
suppa diffSplice \
  -m empirical \
  -p "$PSI_DIR/time_0h.psi" "$PSI_DIR/time_16h.psi" \
  -e "$TPM_DIR/time_0h.tpm" "$TPM_DIR/time_16h.tpm" \
  -i "$FILTERED_DIR/filtered_events_all_replicates.ioe" \
  -gc \
  -pa \
  -a 1000 \
  -o "$DIFF_DIR/0vs16"

# 0h vs 20h
suppa diffSplice \
  -m empirical \
  -p "$PSI_DIR/time_0h.psi" "$PSI_DIR/time_20h.psi" \
  -e "$TPM_DIR/time_0h.tpm" "$TPM_DIR/time_20h.tpm" \
  -i "$FILTERED_DIR/filtered_events_all_replicates.ioe" \
  -gc \
  -pa \
  -a 1000 \
  -o "$DIFF_DIR/0vs20"

# 0h vs 24h
suppa diffSplice \
  -m empirical \
  -p "$PSI_DIR/time_0h.psi" "$PSI_DIR/time_24h.psi" \
  -e "$TPM_DIR/time_0h.tpm" "$TPM_DIR/time_24h.tpm" \
  -i "$FILTERED_DIR/filtered_events_all_replicates.ioe" \
  -gc \
  -pa \
  -a 1000 \
  -o "$DIFF_DIR/0vs24"
```

---

## Helper scripts (repository)
- `scripts/check_tpm_values.py` — validates TPM format and numeric expression values
- `scripts/analyze_all_replicates.py` — filters events keeping only transcripts shared across all replicates
- `scripts/merge_replicates.py` — merges PSI/TPM by timepoint for differential analyses
