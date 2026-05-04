"""
Filter and split sample IDs for train/test sets from a BCF file.

Usage:
    python 01_1_select_samples.py \
        <lab_info.csv> \
        <id_list.txt> \
        <phased.bcf> \
        [output_prefix]

Logic:
    - Extract valid sample IDs from the provided phased BCF file using bcftools.
    - Read the phenotype/ID list and cross-reference with the valid BCF IDs.
    - Retain chronological ordering from the CSV input.
    - Split the overlapping IDs into an 80% training set and a 20% testing set.

Outputs:
    <output_prefix>_train_80.txt  - Training set IDs (plain text)
    <output_prefix>_train_80.indv - Training set IDs (PLINK format)
    <output_prefix>_test_20.txt   - Testing set IDs (plain text)
    <output_prefix>_test_20.indv  - Testing set IDs (PLINK format)
"""

import pandas as pd
import subprocess
import sys
import os

def get_valid_ids_from_bcf(bcf_path):
    """Run bcftools query -l to extract sample IDs from BCF file."""
    result = subprocess.run(
        ["bcftools", "query", "-l", bcf_path],
        capture_output=True, text=True, check=True
    )
    ids = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    print(f"BCF sample IDs extracted: {len(ids)}")
    return set(ids)

def write_outputs(ids, prefix, label):
    """Write both plain txt and plink indvlist format."""
    txt_path   = f"{prefix}_{label}.txt"
    plink_path = f"{prefix}_{label}.indv"

    with open(txt_path, "w") as f:
        f.write("\n".join(ids) + "\n")

    with open(plink_path, "w") as f:
        f.write("\n".join(f"{id} {id}" for id in ids) + "\n")

    print(f"  {txt_path}")
    print(f"  {plink_path}")

def split_twb2_ids(csv_path, id_list_path, bcf_path, output_prefix="split"):
    # Step 1: Get valid IDs from BCF
    valid_bcf_ids = get_valid_ids_from_bcf(bcf_path)

    # Step 2: Load the chronologically sorted CSV
    df = pd.read_csv(csv_path)

    ordered_twb2 = (
        df["TWB2_ID"]
        .dropna()
        .replace("", pd.NA)
        .dropna()
        .tolist()
    )

    # Step 3: Load your input TWB2 ID list
    with open(id_list_path) as f:
        input_ids = [line.strip() for line in f if line.strip()]
    input_set = set(input_ids)

    # Step 4: Filter to IDs that are in BOTH the input list AND the BCF
    valid_input_ids = input_set & valid_bcf_ids
    print(f"Input IDs: {len(input_ids)}, valid after BCF filter: {len(valid_input_ids)}")

    # Step 5: Filter ordered_twb2 to only those present in valid_input_ids (preserves chronological order)
    ordered_input = [twb2_id for twb2_id in ordered_twb2 if twb2_id in valid_input_ids]

    # Warn about IDs not found in CSV or BCF
    found_set = set(ordered_input)
    missing_bcf = [id for id in input_ids if id not in valid_bcf_ids]
    missing_csv = [id for id in input_ids if id not in found_set and id in valid_bcf_ids]

    if missing_bcf:
        print(f"WARNING: {len(missing_bcf)} IDs excluded (not in BCF):")
        for m in missing_bcf:
            print(f"  {m}")
    if missing_csv:
        print(f"WARNING: {len(missing_csv)} IDs excluded (not found in CSV):")
        for m in missing_csv:
            print(f"  {m}")

    # Step 6: 80/20 split
    n = len(ordered_input)
    split_idx = round(n * 0.8)

    train_ids = ordered_input[:split_idx]
    test_ids  = ordered_input[split_idx:]

    # Write outputs
    print(f"\nTotal IDs matched: {n}")
    print(f"Train (80%): {len(train_ids)}")
    write_outputs(train_ids, output_prefix, "train_80")
    print(f"Test  (20%): {len(test_ids)}")
    write_outputs(test_ids, output_prefix, "test_20")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python split_twb2.py <lab_info.csv> <id_list.txt> <phased.bcf> [output_prefix]")
        sys.exit(1)

    csv_path     = sys.argv[1]
    id_list_path = sys.argv[2]
    bcf_path     = sys.argv[3]
    prefix       = sys.argv[4] if len(sys.argv) > 4 else "split"

    split_twb2_ids(csv_path, id_list_path, bcf_path, prefix)
    