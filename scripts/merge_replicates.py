import argparse
from datetime import datetime, timezone
import os

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge PSI/TPM replicate groups by timepoint to enable differential splicing comparisons."
    )

    parser.add_argument("--psi-a", required=True, help="Path to A_psi.psi")
    parser.add_argument("--psi-b", required=True, help="Path to B_psi.psi")
    parser.add_argument("--psi-c", required=True, help="Path to C_psi.psi")

    parser.add_argument("--tpm-a", required=True, help="Path to A_all.tpm")
    parser.add_argument("--tpm-b", required=True, help="Path to B_all.tpm")
    parser.add_argument("--tpm-c", required=True, help="Path to C_all.tpm")

    parser.add_argument("--out-psi-dir", required=True, help="Output directory for time_*.psi files")
    parser.add_argument("--out-tpm-dir", required=True, help="Output directory for time_*.tpm files")

    parser.add_argument(
        "--timepoints",
        default="0h,16h,20h,24h",
        help="Comma-separated timepoint labels matching columns 1..N (default: 0h,16h,20h,24h).",
    )
    parser.add_argument(
        "--replicate-prefixes",
        default="A,B,C",
        help="Comma-separated replicate prefixes for output columns (default: A,B,C).",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()
    print(f"Started: {datetime.now(timezone.utc).isoformat()}")

    os.makedirs(args.out_psi_dir, exist_ok=True)
    os.makedirs(args.out_tpm_dir, exist_ok=True)

    timepoints = [t.strip() for t in args.timepoints.split(",") if t.strip()]
    reps = [r.strip() for r in args.replicate_prefixes.split(",") if r.strip()]
    if reps != ["A", "B", "C"]:
        print("Note: this script currently expects three replicate groups (A, B, C).")

    # Read PSI with index as event_id
    psi_a = pd.read_csv(args.psi_a, sep="\t", index_col=0)
    psi_b = pd.read_csv(args.psi_b, sep="\t", index_col=0)
    psi_c = pd.read_csv(args.psi_c, sep="\t", index_col=0)

    # Read TPM with index as transcript_id
    tpm_a = pd.read_csv(args.tpm_a, sep="\t", index_col=0)
    tpm_b = pd.read_csv(args.tpm_b, sep="\t", index_col=0)
    tpm_c = pd.read_csv(args.tpm_c, sep="\t", index_col=0)

    # Validate number of columns
    n = len(timepoints)
    if psi_a.shape[1] < n or tpm_a.shape[1] < n:
        raise SystemExit("Not enough columns in PSI/TPM inputs for the provided timepoints list.")

    for i, tp in enumerate(timepoints):
        col_idx = i  # 0-based
        # PSI
        out_psi = pd.DataFrame(
            {
                "event_id": psi_a.index,
                f"A_{tp}": psi_a.iloc[:, col_idx],
                f"B_{tp}": psi_b.iloc[:, col_idx],
                f"C_{tp}": psi_c.iloc[:, col_idx],
            }
        )
        out_psi_path = os.path.join(args.out_psi_dir, f"time_{tp}.psi")
        out_psi.to_csv(out_psi_path, sep="\t", index=False)

        # TPM
        out_tpm = pd.DataFrame(
            {
                "transcript_id": tpm_a.index,
                f"A_{tp}": tpm_a.iloc[:, col_idx],
                f"B_{tp}": tpm_b.iloc[:, col_idx],
                f"C_{tp}": tpm_c.iloc[:, col_idx],
            }
        )
        out_tpm_path = os.path.join(args.out_tpm_dir, f"time_{tp}.tpm")
        out_tpm.to_csv(out_tpm_path, sep="\t", index=False)

        print(f"Saved: {out_psi_path}")
        print(f"Saved: {out_tpm_path}")

    print(f"Completed: {datetime.now(timezone.utc).isoformat()}")


if __name__ == "__main__":
    main()
