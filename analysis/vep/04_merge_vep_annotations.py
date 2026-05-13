#!/usr/bin/env python3
"""
Script: 04_merge_vep_annotations.py
Description:
    Merges VEP annotation files (from step 2 and 3) with significant variant
    association tables (from step 1).

    For each combination of phenotype type (quan/bi) and variant type (snv/sv),
    the script loads pre-computed VEP annotations across chr1–22 and left-joins
    them onto the association table by variant ID.

================================================================================
INPUTS:
- VEP annotation TSVs in : ${OUTPUT_DIR}/vep/{snv|sv}/chr{N}.rsq0.8_maf0.005.tsv
- Association tables in  : ${OUTPUT_DIR}/summary_associations/tables/

================================================================================
OUTPUT:
- Annotated TSVs in      : ${OUTPUT_DIR}/summary_associations/tables/
  - quan_sv_vep_annotated.tsv
  - bi_sv_vep_annotated.tsv
  - quan_snv_vep_annotated.tsv
  - bi_snv_vep_annotated.tsv
"""

import os
import pandas as pd
from pathlib import Path

# --------------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------------
base_out_dir = Path(os.environ.get("OUTPUT_DIR", "/staging/biology/u4432941/sv/prs/outputs"))
vep_base     = base_out_dir / "vep"
table_dir    = base_out_dir / "summary_associations" / "tables"

# --------------------------------------------------------------------------------
# Helper functions
# --------------------------------------------------------------------------------
def load_vep_annotations(var_type: str):
    """
    Load and concatenate VEP annotation files across chr1–22 for a given variant type.

    Args:
        var_type : Variant type subdirectory, either "snv" or "sv".

    Returns:
        A single concatenated DataFrame of all available chromosomes.
        Missing chromosome files are skipped with a warning.
    """
    vep_dir = vep_base / var_type
    all_dfs = []

    for chr_num in range(1, 23):
        file_path = vep_dir / f"chr{chr_num}.rsq0.8_maf0.005.tsv"
        if os.path.exists(file_path):
            all_dfs.append(pd.read_csv(file_path, sep="\t", low_memory=False))
        else:
            print(f"[Warning] Missing VEP file: {file_path}")

    if not all_dfs:
        raise FileNotFoundError(f"No VEP annotation files found for var_type='{var_type}' in {vep_dir}")

    return pd.concat(all_dfs, ignore_index=True)

def merge_and_save(vep_df, assoc_file, output_file):
    """
    Left-join a VEP annotation DataFrame onto a variant association table and save the result.

    Variants in the association table without a VEP match are retained (left join),
    with annotation columns filled as NaN.

    Args:
        vep_df      : Combined VEP annotation DataFrame (from load_vep_annotations).
        assoc_file  : Path to the input association TSV (columns must include "ID").
        output_file : Path to write the merged output TSV.
    """
    if not assoc_file.exists():
        print(f"[Warning] Association file not found, skipping: {assoc_file}")
        return

    assoc_df = pd.read_csv(assoc_file, sep="\t")
    merged_df = assoc_df.merge(vep_df, on="ID", how="left")
    merged_df.to_csv(output_file, sep="\t", index=False)
    print(f"[Saved] {output_file}  ({len(merged_df)} rows)")

if __name__ == "__main__":

    # Define all (phenotype_prefix, variant_type, assoc_suffix) combinations to process
    tasks = [
        ("quan", "sv",  "quan_sv_associations.tsv",        "quan_sv_vep_annotated.tsv"),
        ("bi",   "sv",  "bi_sv_associations.tsv",          "bi_sv_vep_annotated.tsv"),
        ("quan", "snv", "quan_snv_associations.tsv",       "quan_snv_vep_annotated.tsv"),
        ("bi",   "snv", "bi_snv_associations.tsv",         "bi_snv_vep_annotated.tsv"),
    ]

    # Load VEP annotations once per variant type to avoid redundant I/O
    print("Loading VEP annotations...")
    vep_cache = {}
    for _, var_type, _, _ in tasks:
        if var_type not in vep_cache:
            vep_cache[var_type] = load_vep_annotations(var_type)
    print(f"Loaded VEP annotation sets: {list(vep_cache.keys())}\n")

    # Merge and save each combination
    print("Merging association files with VEP annotations...")
    for _, var_type, assoc_filename, output_filename in tasks:
        merge_and_save(
            vep_df      = vep_cache[var_type],
            assoc_file  = table_dir / assoc_filename,
            output_file = table_dir / output_filename,
        )

    print("\n--- All tasks completed successfully ---")
