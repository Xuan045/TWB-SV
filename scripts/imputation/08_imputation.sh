#!/usr/bin/bash
#SBATCH -A MST109178
#SBATCH -J imputation
#SBATCH -p ngs92G
#SBATCH -c 14
#SBATCH --mem=92g
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --array=1-22

# ---- Source Configuration ----
source "$(dirname "$0")/../config.sh"

# ---- Environment ----
# Environment is initialized in config.sh
activate_env assoc_env

# ---- Parameters ----
batch='1'
TARGET="sv_snv"
CHR=${SLURM_ARRAY_TASK_ID}

MINIMAC4="${MINIMAC4_BIN}"
PANEL="${IMP_PANEL_PREFIX}/${TARGET}_panel/chr${CHR}.msav"

impdir=${OUT_DIR}/imputation_v${batch}
phasingdir=${impdir}/phasing
outdir=${WORK_DIR:-/work/u4432941/sv}/imputation_v${batch}/${TARGET}_imp
VCF=${phasingdir}/phased.chr${CHR}.bcf

mkdir -p ${outdir}

# ---- Redirect LOG file ----
logdir=${LOG_DIR}/imputation
mkdir -p ${logdir}
logfile=${logdir}/08_imputation_v${batch}_${TARGET}_chr${CHR}.log
exec > "$logfile" 2>&1

# Recode the chr contig name from "1" to "chr1" and etc.
if [ ! -f "$outdir/chr_map.txt" ]; then
  for i in {1..22}; do echo -e "$i\tchr$i"; done > "$outdir/chr_map.txt"
fi
bcftools annotate --rename-chrs "$outdir/chr_map.txt" -o $outdir/chr$CHR.rename.bcf -Ob $VCF
bcftools index $outdir/chr$CHR.rename.bcf

# Run imputation
$MINIMAC4 \
  --threads $SLURM_CPUS_PER_TASK \
  --format GT \
  --all-typed-sites \
  --output $outdir/chr$CHR.dose.vcf.gz \
  --empirical-output $outdir/chr$CHR.empiricalDose.vcf.gz \
  --temp-prefix $outdir \
  $PANEL \
  $outdir/chr$CHR.rename.bcf

# Index the output VCF
bcftools index $outdir/chr$CHR.dose.vcf.gz
bcftools index $outdir/chr$CHR.empiricalDose.vcf.gz

rm $outdir/chr$CHR.rename.bcf*
