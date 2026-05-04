#!/usr/bin/bash
#SBATCH -A MST109178        # Account name/project number
#SBATCH -J POST_IMP_QC      # Job name
#SBATCH -p ngs186G           # Partition Name 
#SBATCH -c 28               
#SBATCH --mem=175G           # memory used
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --array=1-22

# ---- Source Configuration ----
source "$(dirname "$0")/../config.sh"

# ---- Environment ----
# Environment is initialized in config.sh
activate_env assoc_env

# ---- User Configuration ----
batch='2'
TARGET="sv_snv"                  # Imputation panel: "sv_snv" or "snv"
CHR=${SLURM_ARRAY_TASK_ID}

# Step 1 thresholds (broad filtering → PLINK conversion)
STEP1_R2=0.8
STEP1_MAF=0.005
 
# Step 2 thresholds (strict SNV-only filtering → LD pruning varlist)
STEP2_R2=0.9
STEP2_MAF=0.05
 
THREADS=${SLURM_CPUS_PER_TASK}

# Paths
VCF=${WORK_DIR:-/work/u4432941/sv}/imputation_v${batch}/${TARGET}_imp/chr${CHR}.dose.vcf.gz
POST_IMP_DIR="${WORK_DIR:-/work/u4432941/sv}/${TARGET}_post_imp_v${batch}"
PFILE="${POST_IMP_DIR}/chr${CHR}.rsq${STEP1_R2}_maf${STEP1_MAF}"
LD_OUT="${POST_IMP_DIR}/${TARGET}.rsq${STEP2_R2}_maf${STEP2_MAF}_snv"

mkdir -p $POST_IMP_DIR $LD_OUT

# ---- Redirect LOG file ----
logdir=${LOG_DIR}/post_imp_qc
mkdir -p ${logdir}
logfile=${logdir}/09_post_imputation_v${batch}_${TARGET}_chr${CHR}.log
exec > "$logfile" 2>&1

echo "=============================="
echo "Job Name   : $SLURM_JOB_NAME"
echo "Job ID     : $SLURM_JOB_ID"
echo "Node       : $SLURMD_NODENAME"
echo "Start Time : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Working Dir: $(pwd)"
echo "=============================="

# ── Step 1: Broad filtering + PLINK conversion ───
echo "[$(date '+%H:%M:%S')] Step 1 — Filtering VCF (R2>${STEP1_R2}, MAF>${STEP1_MAF}) for chr${CHR}..."
 
bcftools filter \
    -i "R2>${STEP1_R2} && MAF>${STEP1_MAF}" \
    --threads "${THREADS}" \
    -o "${PFILE}.vcf.gz" -Oz \
    "$VCF"
bcftools index --threads "${THREADS}" -f "${PFILE}.vcf.gz"
 
bcftools query -f '%ID\n' "${PFILE}.vcf.gz" \
    > "${PFILE}.varlist"
 
plink2 \
    --vcf "${PFILE}.vcf.gz" \
    --double-id \
    --make-pgen \
    --threads "${THREADS}" \
    --out "$PFILE"
 
echo "[$(date '+%H:%M:%S')] Step 1 done — PLINK pgen written to ${PFILE}"
 
# ── Step 2: Strict SNV filtering + LD pruning ────
echo "[$(date '+%H:%M:%S')] Step 2 — Extracting SNVs (R2>${STEP2_R2}, MAF>${STEP2_MAF}) for chr${CHR}..."
 
VARLIST="${LD_OUT}/chr${CHR}.rsq${STEP2_R2}_maf${STEP2_MAF}_snv.varlist"
LD_PREFIX="${LD_OUT}/chr${CHR}.rsq${STEP2_R2}_maf${STEP2_MAF}_snv.ldpr.geno0.05"
 
bcftools view -v snps --threads "${THREADS}" "$VCF" | \
    bcftools filter -i "R2>${STEP2_R2} && MAF>${STEP2_MAF}" | \
    bcftools query -f '%ID\n' > "$VARLIST"
 
plink2 \
    --pfile "$PFILE" \
    --extract "$VARLIST" \
    --indep-pairwise 500kb 0.9 \
    --geno 0.05 \
    --threads "${THREADS}" \
    --out "$LD_PREFIX"
 
echo "[$(date '+%H:%M:%S')] Step 2 done — LD pruning output at ${LD_PREFIX}"
 
EXIT_CODE=$?
echo "=============================="
echo "End Time   : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Exit Code  : $EXIT_CODE"
echo "=============================="
exit $EXIT_CODE