#!/usr/bin/bash
#SBATCH -A MST109178
#SBATCH -J liftover_hg38
#SBATCH -p ngs92G
#SBATCH -c 14
#SBATCH --mem=92G
#SBATCH -o /dev/null
#SBATCH -e /dev/null

# ---- Redirect LOG file ----
logfile=${LOG_DIR}/00_twb1_liftover.log
exec > "$logfile" 2>&1

# ---- Source Configuration ----
source "$(dirname "$0")/../config.sh"

# ---- Environment ----
# Environment is initialized in config.sh
activate_env vcf_liftover

twb1_hg19=/staging/reserve/jacobhsu/TWB/TWB/microarray/TWB1/27719/Genotyped/TWB1.hg19
chain_file=$CHAIN_FILE
ref_hg38_fa=$FASTA_REF
OUTDIR=${OUT_DIR}/twb1_hg38
mkdir -p $OUTDIR

# Recode to VCF file
echo "Converting PLINK to VCF..."
plink2 --bfile $twb1_hg19 \
    --export vcf id-paste=iid \
    --out ${OUTDIR}/twb1_hg19

Liftover
echo "Running CrossMap VCF LiftOver..."
CrossMap vcf \
    $chain_file \
    ${OUTDIR}/twb1_hg19.vcf \
    $ref_hg38_fa \
    ${OUTDIR}/twb1_hg38.vcf

echo "Sorting and cleaning VCF..."
bcftools sort -m 90G -T ${OUTDIR}/tmp ${OUTDIR}/twb1_hg38.vcf | \
    bcftools view -v snps -m2 -M2 | \
    bcftools filter -e 'ALT=="*"' -o ${OUTDIR}/twb1_hg38_clean.vcf

# Convert VCF back to PLINK format
echo "Converting VCF back to sorted pgen..."
plink2 --vcf ${OUTDIR}/twb1_hg38_clean.vcf \
    --fam ${twb1_hg19}.fam \
    --chr 1-24 \
    --double-id \
    --split-par hg38 \
    --output-chr 26 \
    --sort-vars \
    --make-pgen \
    --out ${OUTDIR}/twb1_hg38

echo "Converting sorted pgen to bed..."
plink2 --pfile ${OUTDIR}/twb1_hg38 \
    --make-bed \
    --out ${OUTDIR}/twb1_hg38

# Remove intermediate VCF files
rm ${OUTDIR}/twb1_hg19.vcf ${OUTDIR}/twb1_hg38.vcf ${OUTDIR}/twb1_hg38_clean.vcf
