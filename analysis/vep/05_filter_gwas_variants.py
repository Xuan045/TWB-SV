#!/usr/bin/env python3
"""
Script: 05_filter_gwas_variants.py
Description:
    Filters VEP-annotated SNV files (chr1–22) to retain only variants with at
    least one non-empty GWAS annotation field, then merges them into a single TSV.

================================================================================
INPUTS:
- VEP annotation TSVs in : ${OUTPUT_DIR}/vep/snv/chr{N}.rsq0.8_maf0.005.tsv

================================================================================
OUTPUT:
- ${OUTPUT_DIR}/summary_associations/tables/snv_gwas_variants.tsv
"""

import os
import pandas as pd
from pathlib import Path

# --------------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------------
base_out_dir = Path(os.environ.get("OUTPUT_DIR", "/staging/biology/u4432941/sv/prs/outputs"))
vep_snv_dir  = base_out_dir / "vep" / "snv"
output_file  = base_out_dir / "summary_associations" / "tables" / "snv_gwas_variants.tsv"

# GWAS-related fields to check
GWAS_FIELDS = [
    "GWAS_accessions", "GWAS_associated_gene", "GWAS_beta_coef",
    "GWAS_odds_ratio", "GWAS_p_value", "GWAS_pmid",
    "GWAS_risk_allele", "GWAS_study"
]

# --------------------------------------------------------------------------------
# Helper functions
# --------------------------------------------------------------------------------
def load_and_filter_gwas(chr_num: int) -> pd.DataFrame | None:
    """
    Load a single chromosome's VEP annotation file and retain only rows
    where at least one GWAS field is not '.'.

    Args:
        chr_num : Chromosome number (1–22).

    Returns:
        Filtered DataFrame, or None if the file is missing or no rows pass the filter.
    """
    file_path = vep_snv_dir / f"chr{chr_num}.rsq0.8_maf0.005.tsv"

    if not file_path.exists():
        print(f"[Warning] Missing file: {file_path}")
        return None

    df   = pd.read_csv(file_path, sep="\t", dtype=str, low_memory=False)
    mask = df[GWAS_FIELDS].apply(lambda row: any(val != "." for val in row), axis=1)

    filtered = df[mask]
    if filtered.empty:
        return None

    return filtered

# --------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------
if __name__ == "__main__":
    filtered_dfs = []

    for chr_num in range(1, 23):
        result = load_and_filter_gwas(chr_num)
        if result is not None:
            filtered_dfs.append(result)

    if not filtered_dfs:
        print("[Error] No GWAS-annotated variants found across all chromosomes.")
    else:
        merged_df = pd.concat(filtered_dfs, ignore_index=True)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        merged_df.to_csv(output_file, sep="\t", index=False)
        print(f"[Saved] {output_file}  ({len(merged_df)} rows)")
