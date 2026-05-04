#!/usr/bin/bash
#SBATCH -A MST109178
#SBATCH -J phasing
#SBATCH -p ngs92G
#SBATCH -c 14
#SBATCH --mem=92g
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --array=1-22

# ---- Source Configuration ----
source "$(dirname "$0")/../config.sh"

# ---- Parameters ----
batch='2'
CHR=${SLURM_ARRAY_TASK_ID}

# ---- Redirect LOG file ----
logdir=${LOG_DIR}/phasing
mkdir -p ${logdir}

logfile=${logdir}/07_phasing_v${batch}_chr${CHR}.log
exec > "$logfile" 2>&1

impdir=${OUT_DIR}/imputation_v${batch}
phasingdir=${impdir}/phasing
VCF=${impdir}/format_vcf/twb${batch}_EAS.qc.chr${CHR}.vcf.gz

mkdir -p ${phasingdir}

script_dir="${SHAPEIT5_DIR}"
GMAP="${GMAP_PREFIX}${SLURM_ARRAY_TASK_ID}.b38.gmap.gz"

${script_dir}/phase_common_static \
    --input $VCF \
    --thread $SLURM_CPUS_PER_TASK \
	--region $CHR \
    --map $GMAP \
    --output ${phasingdir}/phased.chr${CHR}.bcf
