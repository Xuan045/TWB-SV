# GWAS Pipeline: Pre-Imputation QC, Imputation, Regenie, and PRS

This repository contains a comprehensive, production-ready pipeline for performing Genome-Wide Association Studies (GWAS) and Polygenic Risk Score (PRS) calculations. The pipeline is designed to handle both Structural Variants (SV) and Single Nucleotide Variants (SNV), executing quality control, data imputation, association testing, and risk score modeling.

## Project Structure

```
├── config/                  # Pipeline configuration
│   └── phenotype_metadata.csv # Centralized phenotype listing
├── envs/                    # Conda/Micromamba environment specifications
├── reference/               # Lightweight static reference data (e.g. .chain files, LD exclusion intervals)
├── scripts/                 # Core execution scripts
│   ├── config.example.sh    # Configuration template for server-specific paths
│   ├── pre_imp_qc/          # Step 1: Pre-Imputation Quality Control
│   ├── imputation/          # Step 2: Genotype Phasing & Imputation
│   ├── regenie_assoc/       # Step 3: Association Testing using REGENIE
│   ├── prs/                 # Step 4: Polygenic Risk Score modeling via PRSice-2
│   └── utils/               # Utility scripts (e.g. atgc filtering, variant classification)
└── README.md                # This document
```

## Setup Instructions & Dependencies

To ensure maximum portability, **all** paths to environments, tools, and reference data have been completely abstracted out of the execution scripts. You must configure them centrally before running the pipeline.

### 1. Conda/Micromamba Environments
The pipeline relies on four distinct Conda environments to manage Python and R dependencies cleanly. Recreate them using the `.yml` files in the `envs/` directory:

```bash
# Example using micromamba
micromamba create -f envs/vcf_liftover.yml
micromamba create -f envs/assoc_env.yml
micromamba create -f envs/data_env.yml
micromamba create -f envs/r_env.yml
```

### 2. External Standalone Binaries
Some core genetic analysis software are **not** managed by Conda and must be compiled or downloaded independently on your HPC:
- **SHAPEIT5** (`phase_common_static`): Used for haplotype phasing in Step 2.
- **Minimac4**: Used for genotype imputation in Step 2.
- **PRSice-2**: Used for Polygenic Risk Score modeling in Step 4.

### 3. Pipeline Configuration (`config.sh`)

You must tell the pipeline where to find your environments, standalone binaries, and large reference data (like genetic maps, FASTA files, and imputation panels).

1. Copy the example configuration:
   ```bash
   cp scripts/config.example.sh scripts/config.sh
   ```
2. Edit `scripts/config.sh` and replace the placeholder `/path/to/your/...` paths with the actual absolute paths on your server.

## Execution Workflow

The pipeline is split into four distinct modules that must be executed sequentially.

### Step 1: Pre-Imputation QC (`pre_imp_qc/`)
Performs critical quality control before imputation:
- Liftover to hg38 (TWB1 only)
- Sample Call Rate Check
- Sex Discrepancy Check
- Heterozygosity Check
- Identity-by-Descent (IBD) / Relatedness Check
- Principal Component Analysis (PCA)
- **Phenotype Preparation**: Extracts and cleans the required traits (from survey/lab data) mapped in the phenotype metadata CSV, handling outliers and missing values, outputting the final `.txt` files for REGENIE.

### Step 2: Imputation (`imputation/`)
- Prepares PLINK formats for VCF.
- Phases genotypes using SHAPEIT.
- Imputes missing variants.

### Step 3: Association (`regenie_assoc/`)
Executes the REGENIE Step 1 and Step 2 algorithms for both continuous and binary traits defined in `config/phenotype_metadata.csv`.
- Features automated merging of chromosome-level results.
- Generates Manhattan and QQ plots per trait.
- Generates a holistic summary plot across all clinical categories.

### Step 4: Polygenic Risk Scores (`prs/`)
Executes PRSice-2 on the output of REGENIE.
- Splits quantitative and binary targets.
- Handles script-specific exclusions safely (e.g., removing low sample-size traits).

## Phenotype Management

Adding or removing a phenotype is fully centralized. Edit `scripts/config/phenotype_metadata.csv` to add the trait column name, assign its type (`Quantitative` or `Binary`), and give it a clinical `GROUP` (e.g., Cardiac, Liver, Bone). The Bash and R scripts will dynamically parse this CSV to structure their workflows.
