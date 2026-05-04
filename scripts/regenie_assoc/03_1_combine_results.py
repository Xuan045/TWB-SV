"""
Combine REGENIE step 2 results across chromosomes for a single phenotype.

Usage:
    python 03_1_combine_results.py --phenotype <phenotype> --panel <panel>

Logic:
    - Parse all .regenie output files for the specified phenotype and panel.
    - Concatenate the files across all 22 chromosomes.
    - Classify variants using sv_type_count module.
    - Save combined output to OUT_DIR and move original chr-level files to WORK_DIR.

Outputs:
    <OUT_DIR>/<panel>.rsq0.8_maf0.005.regenie_results_prs/<phenotype>.regenie
    <WORK_DIR>/regenie/<panel>/chr<N>_<phenotype>.regenie (archived original files)
"""

import os
import re
import shutil
import argparse
import pandas as pd
from glob import glob
from collections import defaultdict
from importlib.util import spec_from_file_location, module_from_spec

# ── Argument Parsing ──────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--phenotype", type=str, required=True, help="Phenotype to process (one at a time)")
parser.add_argument("--panel", type=str, required=True, choices=["sv_snv", "snv"], help="Panel to process")
args = parser.parse_args()
target_pheno = args.phenotype
target_panel = args.panel

# ── Load helper function to classify variant types ────────────────────────────
project_dir = os.environ.get("PROJECT_DIR", "/staging/biology/u4432941/sv/prs")
sv_type_script = os.path.join(project_dir, "scripts/analysis/02_sv_type_count.py")
spec = spec_from_file_location("sv_type_count", sv_type_script)
sv_type_count = module_from_spec(spec)
spec.loader.exec_module(sv_type_count)
classify_variant = sv_type_count.VariantAnalyzer.classify_variant

# ── Paths ─────────────────────────────────────────────────────────────────────
prefix        = "rsq0.8_maf0.005"
staging_dir   = os.environ.get("OUT_DIR", "/staging/biology/u4432941/sv/prs/outputs")
work_base     = os.environ.get("WORK_DIR", "/work/u4432941/sv")
work_dir      = os.path.join(work_base, "regenie", target_panel)
regenie_dir   = os.path.join(staging_dir, f"{target_panel}_regenie_prs")
output_dir    = os.path.join(staging_dir, f"{target_panel}.{prefix}.regenie_results_prs")
os.makedirs(output_dir, exist_ok=True)
os.makedirs(work_dir, exist_ok=True)

# ── Find and group .regenie files by phenotype ────────────────────────────────
# Expected filename: rsq0.8_maf0.005.regenie_step2_{quan|bi}_chr{N}_{PHENOTYPE}.regenie
pattern   = re.compile(rf"{re.escape(prefix)}\.regenie_step2_(quan|bi)_chr(\d+)_(.+)\.regenie$")
all_files = glob(os.path.join(regenie_dir, "*.regenie"))

phenotype_files = defaultdict(list)
for f in all_files:
    match = pattern.search(os.path.basename(f))
    if match:
        phenotype = match.group(3)
        phenotype_files[phenotype].append(f)

if target_pheno not in phenotype_files:
    raise ValueError(f"Phenotype '{target_pheno}' not found in files under {regenie_dir}")

files = phenotype_files[target_pheno]
print(f"\nMerging {len(files)} files for phenotype: {target_pheno}")
for f in sorted(files):
    print(f"  - {f}")

# ── Combine chr files ─────────────────────────────────────────────────────────
def get_chr_number(filepath):
    match = re.search(r"chr(\d+)", filepath)
    return int(match.group(1)) if match else 999

col_dtypes = {
    'CHROM': str, 'GENPOS': int, 'ID': str, 'ALLELE0': str,
    'ALLELE1': str, 'A1FREQ': float, 'N': int, 'TEST': str,
    'BETA': float, 'SE': float, 'CHISQ': float, 'LOG10P': float, 'EXTRA': str,
}

dfs = []
for f in sorted(files, key=get_chr_number):
    dfs.append(pd.read_csv(f, sep=r"\s+", dtype=col_dtypes))
combined_df = pd.concat(dfs, ignore_index=True)

# ── Parse ID → CHROM, POS, REF, ALT ──────────────────────────────────────────
id_pattern = re.compile(r"^(chr)?(?P<chrom>[\w]+):(?P<pos>\d+):(?P<ref>[^:]+):(?P<alt>.+)$")

def parse_id(id_str):
    match = id_pattern.match(id_str)
    if match:
        return pd.Series({'CHROM': match.group('chrom'), 'POS': int(match.group('pos')),
                          'REF': match.group('ref'),     'ALT': match.group('alt')})
    return pd.Series({'CHROM': None, 'POS': None, 'REF': None, 'ALT': None})

combined_df = pd.concat([combined_df, combined_df['ID'].apply(parse_id)], axis=1)
combined_df["VAR_TYPE"] = combined_df.apply(classify_variant, axis=1)

# ── Save combined result ──────────────────────────────────────────────────────
output_path = os.path.join(output_dir, f"{target_pheno}.regenie")
combined_df.to_csv(output_path, sep="\t", index=False)
print(f"\nSaved combined result: {output_path}")

# ── Verify combined file, then copy chr files to /work and delete originals ───
if not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
    raise RuntimeError("Combined output file is missing or empty — skipping file archival")

n_copied, n_deleted = 0, 0
for f in files:
    dest = os.path.join(work_dir, os.path.basename(f))
    try:
        shutil.copy2(f, dest)
        # Verify copy before deleting
        if os.path.getsize(dest) == os.path.getsize(f):
            os.remove(f)
            n_deleted += 1
        else:
            print(f"Warning: size mismatch after copy, keeping original: {f}")
        n_copied += 1
    except Exception as e:
        print(f"Warning: failed to archive {f}: {e}")

print(f"Archived {n_copied} files to {work_dir}, deleted {n_deleted} originals")
