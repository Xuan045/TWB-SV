#!/bin/bash
#SBATCH -A MST109178
#SBATCH -J check_het
#SBATCH -p ngs92G
#SBATCH -c 14
#SBATCH --mem=92G
#SBATCH -o /dev/null
#SBATCH -e /dev/null

# ---- Source Configuration ----
if [ -n "$SLURM_SUBMIT_DIR" ]; then
    if [ -f "$SLURM_SUBMIT_DIR/../config.sh" ]; then
        source "$SLURM_SUBMIT_DIR/../config.sh"
    elif [ -f "$SLURM_SUBMIT_DIR/scripts/config.sh" ]; then
        source "$SLURM_SUBMIT_DIR/scripts/config.sh"
    fi
else
    source "$(dirname "$0")/../config.sh"
fi

# ---- Environment ----
# Environment is initialized in config.sh
activate_env assoc_env

# ---- Parameters ----
GENO_INIT=0.05
GENO_FINAL=0.02
MIND_FINAL=0.02
batch='2'

# ---- Redirect LOG file ----
logfile=${LOG_DIR}/03_check_het_v${batch}.log
exec > "$logfile" 2>&1

preqcdir=${OUT_DIR}/pre_imp_qc_v${batch}
BFILE=${preqcdir}/twb${batch}_geno${GENO_INIT}_mind${MIND_FINAL}_geno${GENO_FINAL}_dedup
PARA=${preqcdir}/twb${batch}

# Perform LD pruning first
plink2 \
    --bfile ${BFILE} \
    --autosome \
    --maf 0.05 \
    --snps-only just-acgt \
    --indep-pairwise 200 100 0.1 \
    --out ${PARA}.ldpr

# Estimate autosomal heterozygosity rate/inbreeding coef
plink \
    --bfile ${BFILE} \
    --extract ${PARA}.ldpr.prune.in \
    --het \
    --out ${PARA}.inbr


# Calculate hetorozygosity rate (i.e., the proportion of heterozygous genotypes for a given individual.)
awk 'NR>1{print ($5-$3)/$5}' ${PARA}.inbr.het | sed '1i HetRate' | \
    paste ${PARA}.inbr.het -> ${PARA}.inbr.hetrate

# Plot
run_in_env r_env Rscript \
    03_check_het.r $preqcdir twb${batch}

# Remove individuals with excessive heterozygosity rate
sort -u ${PARA}_sex_mismatch_F0.20_M0.80.indvlist ${PARA}_het_outlier_6sd.indvlist > ${PARA}_sex_het.indvlist
plink2 \
    --bfile ${BFILE} \
    --remove ${PARA}_sex_het.indvlist \
    --make-bed \
    --out ${PARA}_preqc_sexcheck_hetcheck
    