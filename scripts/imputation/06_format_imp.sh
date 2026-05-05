#!/usr/bin/bash
#SBATCH -A MST109178
#SBATCH -J format_imp
#SBATCH -p ngs92G
#SBATCH -c 14
#SBATCH --mem=92g
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
batch='2'

# ---- Redirect LOG file ----
logfile=${LOG_DIR}/06_format_imp_v${batch}.log
exec > "$logfile" 2>&1

preqcdir=${OUT_DIR}/pre_imp_qc_v${batch}
pcadir=${OUT_DIR}/pca_v${batch}
impdir=${OUT_DIR}/imputation_v${batch}/format_vcf
BFILE=${pcadir}/twb${batch}_preqc_sexcheck_hetcheck_EAS
PARA=${impdir}/twb${batch}

mkdir -p ${impdir}

plink2 \
    --bfile ${BFILE} \
    --hwe 1e-5 0.001 \
    --make-bed \
    --out ${PARA}_EAS.qc

for chr in {1..22}; do
  plink2 \
    --bfile ${PARA}_EAS.qc \
    --chr $chr \
    --recode vcf id-paste=iid bgz \
    --out ${PARA}_EAS.qc.chr${chr}

    # Add AC and AN annotations
    bcftools +fill-tags ${PARA}_EAS.qc.chr${chr}.vcf.gz -- -t AC,AN | \
        bgzip -c > temp && mv temp ${PARA}_EAS.qc.chr${chr}.vcf.gz

    bcftools index ${PARA}_EAS.qc.chr${chr}.vcf.gz
done
