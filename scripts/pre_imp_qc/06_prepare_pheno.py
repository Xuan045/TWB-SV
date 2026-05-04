"""
Usage:
  python3 06_prepare_pheno.py

Logic:
  1. Reads configuration paths (SURVEY_DATA, LAB_DATA, PHENO_META) from environment.
  2. Dynamically loads the official phenotype list (Quantitative & Binary) from PHENO_META.
  3. Maps the official trait names to raw survey column names (e.g., adding '_SELF' if it exists).
  4. Merges survey and lab test data on TWB1_ID and TWB2_ID.
  5. Cleans data (removes outliers, replaces -9 with NA).
  6. Joins with PCA eigenvectors from Step 05.
  7. Outputs twb1_pheno_covar.txt and twb2_pheno_covar.txt for REGENIE.

Outputs:
  - ${PROJECT_DIR}/pheno_summary/twb1_pheno_covar.txt
  - ${PROJECT_DIR}/pheno_summary/twb2_pheno_covar.txt
  - Distribution PDFs and Summary TSVs in the same folder.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf_backend

# ── Paths from Environment ───────────────────────────────────────────────────
bi_pheno_path   = os.environ.get("SURVEY_DATA")
quan_pheno_path = os.environ.get("LAB_DATA")
pheno_meta_path = os.environ.get("PHENO_META")
project_dir     = os.environ.get("PROJECT_DIR")
out_dir         = os.environ.get("OUT_DIR")

if not all([bi_pheno_path, quan_pheno_path, pheno_meta_path, project_dir, out_dir]):
    raise ValueError("Missing essential environment variables (SURVEY_DATA, LAB_DATA, PHENO_META, PROJECT_DIR, OUT_DIR). Did you source config.sh?")

twb1_pca_path   = os.path.join(out_dir, "pca_v1", "twb1_pca.eigenvec")
twb2_pca_path   = os.path.join(out_dir, "pca_v2", "twb2_pca.eigenvec")
pheno_outdir    = os.path.join(project_dir, "pheno_summary")
os.makedirs(pheno_outdir, exist_ok=True)

# ── Load Phenotype Metadata ──────────────────────────────────────────────────
print(f"Loading metadata from: {pheno_meta_path}")
meta_df = pd.read_csv(pheno_meta_path)
binary_traits = meta_df[meta_df['TYPE'] == 'Binary']['PHENOTYPE'].tolist()
quant_traits = meta_df[meta_df['TYPE'] == 'Quantitative']['PHENOTYPE'].tolist()

print(f"Loaded {len(binary_traits)} binary and {len(quant_traits)} quantitative traits.")

# ── Dynamic Column Mapping ───────────────────────────────────────────────────
# Some raw columns have different names than our official phenotype names
trait_map = {"FEV1": "FEV10", "FEV1_FVC": "FEV10_FVC"}

bi_header = pd.read_csv(bi_pheno_path, nrows=0).columns.tolist()
quan_header = pd.read_csv(quan_pheno_path, nrows=0).columns.tolist()

bi_raw_cols = ["TWB1_ID", "TWB2_ID", "AGE", "SEX"]
quan_raw_cols = ["TWB1_ID", "TWB2_ID"]

# Map Binary Traits
for t in binary_traits:
    raw_t = trait_map.get(t, t)
    if raw_t + "_SELF" in bi_header:
        bi_raw_cols.append(raw_t + "_SELF")
    elif raw_t in bi_header:
        bi_raw_cols.append(raw_t)
    elif raw_t in quan_header: # Edge case: some binary might be in lab data
        pass
    else:
        print(f"Warning: Binary trait '{t}' (mapped as '{raw_t}') not found in survey headers.")

# Map Quantitative Traits
for t in quant_traits:
    raw_t = trait_map.get(t, t)
    if raw_t + "_SELF" in quan_header:
        quan_raw_cols.append(raw_t + "_SELF")
    elif raw_t in quan_header:
        quan_raw_cols.append(raw_t)
    elif raw_t in bi_header: # Edge case: some quant might be in survey data
        if raw_t + "_SELF" not in bi_raw_cols and raw_t not in bi_raw_cols:
            bi_raw_cols.append(raw_t + "_SELF" if raw_t + "_SELF" in bi_header else raw_t)
    else:
        print(f"Warning: Quantitative trait '{t}' (mapped as '{raw_t}') not found in lab headers.")

pca_cols = ["IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"]

# ── Functions ─────────────────────────────────────────────────────────────────
def load_in_chunks(path, cols, chunksize=10000):
    """Load CSV in chunks, dropping rows where BOTH TWB1_ID and TWB2_ID are missing."""
    chunks = [chunk.dropna(subset=['TWB1_ID', 'TWB2_ID'], how='all')
              for chunk in pd.read_csv(path, usecols=cols, chunksize=chunksize, low_memory=False)]
    return pd.concat(chunks, ignore_index=True)

def clean_df(df, id_col, pca_df, afp_threshold=100):
    """Merge PCA, strip _SELF suffix, convert -9 to NA, recode SEX, remove AFP outliers, winsorize."""
    df = pd.merge(df, pca_df, left_on=id_col, right_on="IID")
    df.columns = [c.replace('_SELF', '') for c in df.columns]
    df = df.replace(-9, pd.NA)

    # Remove AFP values likely due to pregnancy (females above threshold)
    if 'AFP' in df.columns and 'SEX' in df.columns:
        female_high_afp = (df['SEX'] == 2) & (df['AFP'] > afp_threshold)
        print(f"  Females with AFP > {afp_threshold}: {female_high_afp.sum()} set to NA")
        df.loc[female_high_afp, 'AFP'] = pd.NA

    # Winsorize: set values beyond 6 SD to NA
    for trait in quant_traits:
        if trait not in df.columns:
            continue
        mean, std = df[trait].mean(), df[trait].std()
        mask = (df[trait] < mean - 6 * std) | (df[trait] > mean + 6 * std)
        if mask.sum() > 0:
            print(f"  {trait}: {mask.sum()} outliers (>6 SD) set to NA")
        df.loc[mask, trait] = pd.NA

    return df

def save_pheno(df, id_col, path):
    """Save phenotype file with FID and IID columns as required by REGENIE."""
    out = df.copy()
    out = out.drop(columns=['FID', 'IID'], errors='ignore')
    out.insert(0, 'FID', out[id_col])
    out.insert(1, 'IID', out[id_col])
    out.to_csv(path, sep="\t", index=False, na_rep='NA')

def print_summary(df, label):
    """Print sample count and sex breakdown."""
    if 'SEX' in df.columns:
        sex = df['SEX'].value_counts().reindex([1, 2], fill_value=0)
        print(f"\n── {label} ──────────────────────────")
        print(f"  Total   : {len(df)}")
        print(f"  Male (1): {sex[1]}")
        print(f"  Female (2): {sex[2]}")
        print(f"  SEX NA  : {df['SEX'].isna().sum()}")

def binary_summary(df, traits, label):
    """Summarize control/case/NA counts for binary traits."""
    rows = []
    for trait in traits:
        if trait not in df.columns:
            continue
        col = df[trait]
        rows.append({'Phenotype': trait,
                     'Control': (col == 0).sum(),
                     'Case':    (col == 1).sum(),
                     'NA':      col.isna().sum()})
    return pd.DataFrame(rows)

def plot_quant_distributions(df, traits, pdf_path, label):
    """Save histogram for each quantitative trait to a single PDF."""
    with pdf_backend.PdfPages(pdf_path) as pdf:
        for trait in traits:
            if trait not in df.columns:
                continue
            fig, ax = plt.subplots(figsize=(7, 5))
            ax.hist(df[trait].dropna(), bins=50, color="#69b3a2", edgecolor="black")
            ax.set(title=f"[{label}] Distribution of {trait}", xlabel=trait, ylabel="Count")
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

# ── Load raw phenotype data ───────────────────────────────────────────────────
print("\nLoading binary phenotype file...")
bi_pheno_df   = load_in_chunks(bi_pheno_path,   bi_raw_cols)
print("Loading quantitative phenotype file...")
quan_pheno_df = load_in_chunks(quan_pheno_path, quan_raw_cols)

# Merge binary + quantitative on both ID columns
merged_all = pd.merge(bi_pheno_df, quan_pheno_df, on=["TWB1_ID", "TWB2_ID"])
merged_all = merged_all.rename(columns={"FEV10": "FEV1", "FEV10_FVC": "FEV1_FVC"})
print(f"Renamed FEV10 to FEV1. Total merged records: {len(merged_all)}")

# ── Split into TWB1-only and TWB2-only subsets ────────────────────────────────
twb2_raw = merged_all.dropna(subset=['TWB2_ID']).copy()
twb1_raw = merged_all[merged_all['TWB2_ID'].isna()].dropna(subset=['TWB1_ID']).copy()

print(f"TWB2 records: {len(twb2_raw)}")
print(f"TWB1-only records: {len(twb1_raw)}")

# ── Load PCA files ────────────────────────────────────────────────────────────
print("\nLoading PCA eigenvectors...")
twb1_pca = pd.read_csv(twb1_pca_path, sep=r"\s+", usecols=pca_cols)
twb2_pca = pd.read_csv(twb2_pca_path, sep=r"\s+", usecols=pca_cols)

# ── Clean each cohort ─────────────────────────────────────────────────────────
print("\nCleaning TWB1...")
twb1_df = clean_df(twb1_raw, id_col='TWB1_ID', pca_df=twb1_pca)

print("\nCleaning TWB2...")
twb2_df = clean_df(twb2_raw, id_col='TWB2_ID', pca_df=twb2_pca)

# ── Print summaries ───────────────────────────────────────────────────────────
print_summary(twb1_df, "TWB1")
print_summary(twb2_df, "TWB2")

# ── Save individual cohort files ──────────────────────────────────────────────
save_pheno(twb1_df, 'TWB1_ID', os.path.join(pheno_outdir, "twb1_pheno_covar.txt"))
save_pheno(twb2_df, 'TWB2_ID', os.path.join(pheno_outdir, "twb2_pheno_covar.txt"))
print(f"\nSaved REGENIE inputs to: {pheno_outdir}")

# ── Save merged file ──────────────────────────────────────────────────────────
merged_df = pd.concat([twb1_df, twb2_df], ignore_index=True)
merged_df.to_csv(os.path.join(pheno_outdir, "merged_pheno_covar.txt"), sep="\t", index=False)

# ── Binary traits summary ─────────────────────────────────────────────────────
print("\nGenerating binary traits summaries...")
for df, label in [(twb1_df, "twb1"), (twb2_df, "twb2"), (merged_df, "merged")]:
    binary_summary(df, binary_traits, label).to_csv(
        os.path.join(pheno_outdir, f"{label}_binary_traits_summary.tsv"), sep="\t", index=False)

# ── Quantitative traits distribution plots ────────────────────────────────────
print("Generating quantitative traits distribution PDFs...")
for df, label in [(twb1_df, "twb1"), (twb2_df, "twb2"), (merged_df, "merged")]:
    plot_quant_distributions(
        df, quant_traits,
        os.path.join(pheno_outdir, f"{label}_quant_traits_distributions.pdf"),
        label.upper()
    )

print("Phenotype preparation successfully completed!")
