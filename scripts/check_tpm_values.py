import argparse
import glob
import os
import pandas as pd


def check_file_format(filename: str) -> pd.DataFrame | None:
    try:
        df = pd.read_csv(filename, sep="\t", header=None)

        num_columns = len(df.columns)
        if num_columns != 2:
            print(f"WARNING: {filename} has {num_columns} columns (expected 2).")
            print("First rows:")
            print(df.head())
            return None

        df.columns = ["Transcript_ID", "Expression_Value"]

        non_numeric = pd.to_numeric(df["Expression_Value"], errors="coerce").isna()
        if non_numeric.any():
            print(f"Non-numeric values found in {filename}:")
            print(df[non_numeric])
        else:
            print(f"OK: all expression values are numeric in {filename}.")

        return df

    except Exception as e:
        print(f"Error processing {filename}: {e}")
        return None


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate TPM files: 2 columns and numeric expression values.")
    parser.add_argument("--input-dir", required=True, help="Directory containing TPM files.")
    parser.add_argument("--pattern", default="*.tpm", help="Glob pattern (default: *.tpm).")
    args = parser.parse_args()

    input_dir = args.input_dir
    pattern = args.pattern

    files = sorted(glob.glob(os.path.join(input_dir, pattern)))
    if not files:
        raise SystemExit(f"No files found in {input_dir} matching pattern: {pattern}")

    for filename in files:
        print(f"\nAnalyzing {os.path.basename(filename)} ...")
        df = check_file_format(filename)
        if df is not None:
            print(f"Rows: {len(df)}")
            print(f"Expression range: {df['Expression_Value'].min()} - {df['Expression_Value'].max()}")


if __name__ == "__main__":
    main()
