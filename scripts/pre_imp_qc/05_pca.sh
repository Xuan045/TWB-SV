#!/usr/bin/bash
#SBATCH -A MST109178
#SBATCH -J pca
#SBATCH -p ngs372G
#SBATCH -c 56
#SBATCH --mem=350g
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
batch='1'
pred_prob=0.8

# ---- Redirect LOG file ----
logfile=${LOG_DIR}/05_pca_v${batch}.log
exec > "$logfile" 2>&1

preqcdir=${OUT_DIR}/pre_imp_qc_v${batch}
pcadir=${OUT_DIR}/pca_v${batch}
BFILE=${preqcdir}/twb${batch}_preqc_sexcheck_hetcheck
PARA=${pcadir}/twb${batch}
PED=${ONEKGP_PED_TXT}

mkdir -p ${pcadir}

# ── Study-only PCA ────────────────────────────────────────────────────
run_in_env data_env \
    python3 ${PROJECT_DIR}/scripts/utils/find_atgc_snps.py ${BFILE}.bim > ${PARA}_atgc.snplist

awk 'NR==FNR{a[$1];next} !($2 in a) {print $2}' \
    ${PARA}_atgc.snplist ${BFILE}.bim > ${PARA}_nonatgc.snplist

plink \
    --bfile ${BFILE} \
    --autosome \
    --geno 0.02 \
    --maf 0.05 \
    --snps-only just-acgt \
    --extract ${PARA}_nonatgc.snplist \
    --indep-pairwise 200 100 0.1 \
    --out ${PARA}_ldpr

plink2 \
    --bfile ${BFILE} \
    --extract ${PARA}_ldpr.prune.in \
    --pca approx 20 \
    --out ${PARA}_pca

# ── Ref PCA ───────────────────────────────────────────────────────────
ref=${ONEKGP_REF_PREFIX}

awk '{$2 = "chr"$2; print}' ${ref}.bim > ${PARA}_ref_chr.bim
ln -sf ${ref}.bed ${PARA}_ref_chr.bed
ln -sf ${ref}.fam ${PARA}_ref_chr.fam

awk 'NR==FNR{a[$2];next} ($2 in a){print $2}' \
    ${BFILE}.bim ${PARA}_ref_chr.bim > ${PARA}_ref_comm.snplist

plink2 \
    --bfile ${BFILE} \
    --extract ${PARA}_ref_comm.snplist \
    --make-bed \
    --out ${PARA}_comm

plink2 \
    --bfile ${PARA}_ref_chr \
    --extract ${PARA}_ref_comm.snplist \
    --make-bed \
    --out ${PARA}_ref_comm

plink \
    --bfile ${PARA}_comm \
    --keep-allele-order \
    --bmerge ${PARA}_ref_comm \
    --make-bed \
    --out ${PARA}_ref

rm ${PARA}_ref_chr.bim ${PARA}_ref_chr.bed ${PARA}_ref_chr.fam

run_in_env data_env \
    python3 ${PROJECT_DIR}/scripts/utils/find_atgc_snps.py ${PARA}_ref.bim > ${PARA}_ref_atgc.snplist

awk 'NR==FNR{a[$1];next} !($2 in a) {print $2}' \
    ${PARA}_ref_atgc.snplist ${PARA}_ref.bim > ${PARA}_ref_nonatgc.snplist

plink \
    --bfile ${PARA}_ref \
    --autosome \
    --geno 0.02 \
    --maf 0.05 \
    --snps-only just-acgt \
    --extract ${PARA}_ref_nonatgc.snplist \
    --exclude range $LD_EXCLUDE \
    --indep-pairwise 200 100 0.1 \
    --out ${PARA}_ref_ldpr

plink2 \
    --bfile ${PARA}_ref \
    --extract ${PARA}_ref_ldpr.prune.in \
    --pca approx 20 \
    --out ${PARA}_ref_pca

# ── R: study-only + ref PCA plot + RF prediction ───────────────────────────
run_in_env r_env \
    Rscript 05_pca.r ${pcadir} twb${batch} TRUE \
    ${PROJECT_DIR}/scripts/utils/color_code.R

# ── Extract individuals pass PCA check (not outliers and is EAS) ────────────────
plink2 \
    --bfile ${BFILE} \
    --keep ${PARA}_ref_predpop${pred_prob}_EAS.indvlist \
    --make-bed \
    --out ${PARA}_preqc_sexcheck_hetcheck_EAS

# ── EAS-only PCA（after RF prediction）──────────────────────────────────
# Get EAS 1kGP samples from pedigree
awk 'NR>1 && $NF=="EAS" {print $2, $2}' ${PED} \
    > ${PARA}_1kg_EAS.indvlist

# Merge EAS 1kGP + TWB EAS predicted samples
cat ${PARA}_1kg_EAS.indvlist \
    ${PARA}_ref_predpop${pred_prob}_EAS.indvlist \
    > ${PARA}_EAS_keep.indvlist

plink2 \
    --bfile ${PARA}_ref \
    --keep ${PARA}_EAS_keep.indvlist \
    --extract ${PARA}_ref_ldpr.prune.in \
    --pca approx 20 \
    --out ${PARA}_EAS_pca

# ── R: EAS-only plot ──────────────────────────────────────────────────
run_in_env r_env \
    Rscript 05_pca.r ${pcadir} twb${batch} TRUE \
    ${PROJECT_DIR}/scripts/utils/color_code.R
