import argparse
import os
from datetime import datetime, timezone

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter events.ioe keeping only events whose transcripts are present in all replicate TPM matrices."
    )
    parser.add_argument("--events-ioe", required=True, help="Path to combined events.ioe file.")
    parser.add_argument("--tpm-a", required=True, help="Path to A_all.tpm (joined TPM matrix).")
    parser.add_argument("--tpm-b", required=True, help="Path to B_all.tpm (joined TPM matrix).")
    parser.add_argument("--tpm-c", required=True, help="Path to C_all.tpm (joined TPM matrix).")
    parser.add_argument("--out-dir", required=True, help="Output directory.")
    parser.add_argument(
        "--out-name",
        default="filtered_events_all_replicates.ioe",
        help="Output file name (default: filtered_events_all_replicates.ioe).",
    )
    return parser.parse_args()


def check_transcripts(row: pd.Series, available_transcripts: set[str]) -> bool:
    alt_trans = set(str(row["alternative_transcripts"]).split(","))
    total_trans = set(str(row["total_transcripts"]).split(","))
    return all(t in available_transcripts for t in alt_trans.union(total_trans))


def main() -> None:
    args = parse_args()
    print(f"Started: {datetime.now(timezone.utc).isoformat()}")

    os.makedirs(args.out_dir, exist_ok=True)

    # Read events
    events = pd.read_csv(args.events_ioe, sep="\t")
    print(f"Total events in input: {len(events)}")

    # Collect transcripts from events
    event_transcripts: set[str] = set()
    for _, row in events.iterrows():
        event_transcripts.update(str(row["alternative_transcripts"]).split(","))
        event_transcripts.update(str(row["total_transcripts"]).split(","))

    print(f"Total transcripts referenced by events: {len(event_transcripts)}")

    # Read expression matrices and collect transcript IDs (index)
    expression_files = {"A": args.tpm_a, "B": args.tpm_b, "C": args.tpm_c}
    expression_transcripts: dict[str, set[str]] = {}

    for label, path in expression_files.items():
        expr = pd.read_csv(path, sep="\t", index_col=0)
        expression_transcripts[label] = set(expr.index)
        print(f"{label}: transcripts in TPM matrix: {len(expression_transcripts[label])}")

    common_transcripts = set.intersection(*expression_transcripts.values())
    print(f"Transcripts present in all replicates: {len(common_transcripts)}")

    # Filter events
    events["valid"] = events.apply(lambda row: check_transcripts(row, common_transcripts), axis=1)
    valid_events = events[events["valid"]].drop(columns=["valid"])

    out_path = os.path.join(args.out_dir, args.out_name)
    valid_events.to_csv(out_path, sep="\t", index=False)
    print(f"Saved filtered events: {out_path}")
    print(f"Events after filtering: {len(valid_events)}")

    # Optional: event type summary
    def event_type(event_id: str) -> str:
        try:
            return event_id.split(";")[1].split(":")[0]
        except Exception:
            return "NA"

    events_types = events["event_id"].map(event_type).value_counts()
    valid_types = valid_events["event_id"].map(event_type).value_counts()

    print("\nEvent types (original):")
    print(events_types.to_string())
    print("\nEvent types (filtered):")
    print(valid_types.to_string())

    print(f"\nCompleted: {datetime.now(timezone.utc).isoformat()}")


if __name__ == "__main__":
    main()
