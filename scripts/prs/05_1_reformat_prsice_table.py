"""
Reformat REGENIE summary statistics for PRSice-2 input.

Usage:
    python 05_1_reformat_prsice_table.py \
        --sumstats <combined_regenie_file> \
        --out <output_formatted_file>

Logic:
    - Read the combined REGENIE summary statistics file.
    - Check for all required columns (ID, ALLELE1, ALLELE0, BETA, SE, LOG10P).
    - Convert LOG10P values back to raw P-values.
    - Filter out rows with invalid standard errors (SE <= 0) or P-values <= 0.
    - Rename and select columns to match the PRSice-2 standard format.

Outputs:
    <out> - A tab-separated file formatted for PRSice-2 containing: SNP, CHR, POS, A1, A2, BETA, SE, P.
"""

import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Reformat REGENIE sumstats for PRSice-2 input")
    parser.add_argument("--sumstats", required=True, help="Path to combined REGENIE sumstats file")
    parser.add_argument("--out",      required=True, help="Path to output file")
    return parser.parse_args()

def main():
    args = parse_args()

    df = pd.read_csv(args.sumstats, sep=r"\s+")

    required_cols = ["ID", "ALLELE1", "ALLELE0", "BETA", "SE", "LOG10P"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns: {missing}")

    # Convert LOG10P back to P
    df["P"] = 10 ** (-df["LOG10P"])

    # PRSice-2 wants: SNP, A1, A2, BETA, SE, P
    out_df = pd.DataFrame({
        "SNP":  df["ID"],
        "CHR":  df["CHROM"],
        "POS":  df["GENPOS"],
        "A1":   df["ALLELE1"],
        "A2":   df["ALLELE0"],
        "BETA": df["BETA"],
        "SE":   df["SE"],
        "P":    df["P"]
    })

    n_before = len(out_df)
    out_df = out_df.dropna()
    out_df = out_df[out_df["SE"] > 0]
    out_df = out_df[out_df["P"] > 0]
    n_dropped = n_before - len(out_df)
    if n_dropped > 0:
        print(f"Warning: dropped {n_dropped} rows with missing/invalid values")

    out_df.to_csv(args.out, sep="\t", index=False)
    print(f"Done: {len(out_df)} SNPs written to {args.out}")

if __name__ == "__main__":
    main()
    