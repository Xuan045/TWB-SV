#!/usr/bin/env python3
"""
Script: 01_separate_snv_sv.py
Description:
    Summarizes significant variant associations (P < 5e-8) from REGENIE results
    for both quantitative ("quan") and binary ("bi") phenotypes.

INPUTS:
- Phenotype metadata CSV ($PHENO_META)
- Filtered REGENIE files in: ${OUTPUT_DIR}/sv_snv.rsq0.8_maf0.005.regenie_results/filtered_5e-8/

OUTPUT STRUCTURE:
All outputs are organized within: ${OUTPUT_DIR}/summary_associations/

1. /tables/ (Summary TSVs)
   - {quan|bi}_snv_indel_associations.tsv : SNV/INDELs and their associated traits.
   - {quan|bi}_sv_associations.tsv        : SVs (INS/DEL) and their associated traits.
   (Columns: ID, associated_phenotypes)

2. /varlists/ (Chromosome-separated lists for PLINK/bcftools)
   - /quan/chr{N}_{snv|sv}.varlist
   - /bi/chr{N}_{snv|sv}.varlist
   (Contains unique variant IDs, one per line)
"""

import os
import pandas as pd
from pathlib import Path

# Get output directory from environment variable
base_out_dir = Path(os.environ.get("OUTPUT_DIR", "/staging/biology/u4432941/sv/prs/outputs"))
input_dir = base_out_dir / "sv_snv.rsq0.8_maf0.005.regenie_results" / "filtered_5e-8"

summary_dir = base_out_dir / "summary_associations"
table_dir = summary_dir / "tables"
varlist_root = summary_dir / "varlists"

table_dir.mkdir(parents=True, exist_ok=True)
varlist_root.mkdir(parents=True, exist_ok=True)

# ---------------------------------------
# Helper: Process a single phenotype set
# ---------------------------------------
def process_phenotype_set(phenotypes, prefix):
    """
    Processes a list of phenotypes, aggregates variant-trait associations,
    and writes summary tables and variant lists.
    """
    print(f"--- Starting {prefix.upper()} set ---")

    # Create sub-directory for varlists (e.g., summary_associations/varlists/quan/)
    prefix_varlist_dir = varlist_root / prefix
    prefix_varlist_dir.mkdir(parents=True, exist_ok=True)

    all_snv = []
    all_sv  = []

    for pheno in phenotypes:
        file_path = input_dir / f"{pheno}.filtered.regenie"

        if not os.path.exists(file_path):
            print(f"[Warning] File not found: {file_path}")
            continue

        try:
            df = pd.read_csv(file_path, sep="\t", usecols=["ID", "VAR_TYPE"])
        except Exception as e:
            print(f"[Error] {pheno}: {e}")
            continue

        df["VAR_TYPE"] = df["VAR_TYPE"].str.upper()
        df["pheno"] = pheno

        all_snv.append(df[df["VAR_TYPE"].isin({"SNV", "INDEL"})][["ID", "pheno"]])
        all_sv.append(df[df["VAR_TYPE"].isin({"DEL", "INS"})][["ID", "pheno"]])

    # ---------------------------------------
    # Helper: Save summary TSV and varlists
    # ---------------------------------------
    def save_summary(all_dfs, suffix):
        if not all_dfs:
            print(f"[Warning] No {suffix} data found for {prefix} set")
            return

        combined = pd.concat(all_dfs, ignore_index=True)

        # Aggregate: one row per variant, associated phenotypes as comma-separated string
        summary_df = (
            combined.groupby("ID")["pheno"]
            .apply(lambda x: ",".join(sorted(set(x))))
            .reset_index()
            .rename(columns={"pheno": "associated_phenotypes"})
        )

        # Save summary TSV
        tsv_path = table_dir / f"{prefix}_{suffix}_associations.tsv"
        summary_df.to_csv(tsv_path, sep="\t", index=False)
        print(f"[Saved] Table: {tsv_path}")

        # Write per-chromosome varlists
        chr_groups = summary_df["ID"].str.split(":").str[0]
        for chrom, group_df in summary_df.groupby(chr_groups):
            varlist_path = prefix_varlist_dir / f"{chrom}_{suffix}.varlist"
            with open(varlist_path, "w") as f:
                f.write("\n".join(sorted(group_df["ID"].unique())) + "\n")

        print(f"[Saved] {chr_groups.nunique()} varlist files to {prefix_varlist_dir}")

    save_summary(all_snv, "snv")
    save_summary(all_sv, "sv")

# ---------------------------------------
# Run for both quan and bi phenotypes
# ---------------------------------------
if __name__ == "__main__":
    pheno_meta_file = Path(os.environ.get("PHENO_META", "/staging/biology/u4432941/sv/prs/scripts/config/phenotype_metadata.csv"))

    if not pheno_meta_file.exists():
        raise FileNotFoundError(f"Metadata not found: {pheno_meta_file}")

    pheno_df = pd.read_csv(pheno_meta_file)

    quan_list = pheno_df[pheno_df["TYPE"] == "Quantitative"]["PHENOTYPE"].tolist()
    bi_list = pheno_df[pheno_df["TYPE"] == "Binary"]["PHENOTYPE"].tolist()

    if quan_list:
        process_phenotype_set(quan_list, "quan")
    
    if bi_list:
        process_phenotype_set(bi_list, "bi")

    print("\n--- All tasks completed successfully ---")
