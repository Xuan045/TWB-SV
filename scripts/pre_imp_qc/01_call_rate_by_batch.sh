#!/usr/bin/bash
#SBATCH -A MST109178
#SBATCH -J call_rate
#SBATCH -p ngs92G
#SBATCH -c 14
#SBATCH --mem=92G
#SBATCH -o /dev/null
#SBATCH -e /dev/null

# ---- Source Configuration ----
source "$(dirname "$0")/../config.sh"

# ---- Environment ----
# Environment is initialized in config.sh
activate_env assoc_env

# ---- Parameters ----
REF=$FASTA_REF

GENO_INIT=0.05
GENO_FINAL=0.02
MIND_FINAL=0.02
batch='1'

if [ "$batch" == "1" ]; then
    datadir=${OUT_DIR}/twb1_hg38
    prefix=twb1_hg38
else
    datadir=/staging/reserve/jacobhsu/TWB/TWB/microarray/TWB2_genotyped_120163
    prefix=TWB2.hg38.v4
fi

# ---- Redirect LOG file ----
logfile=${LOG_DIR}/01_call_rate_by_batch_v${batch}.log
exec > "$logfile" 2>&1

preqcdir=${OUT_DIR}/pre_imp_qc_v${batch}
mkdir -p $preqcdir $LOG_DIR

# 1. Output a list of SNPs with call rate >0.95
plink \
    --bfile ${datadir}/${prefix} \
    --allow-extra-chr \
    --chr 1-25 \
    --geno $GENO_INIT \
    --missing \
    --write-snplist \
    --out ${preqcdir}/twb${batch}_geno${GENO_INIT}

# 2. Output a list of samples with call rate >=0.98
awk 'NR>1 {if($6<'$MIND_FINAL') print $1,$2}' ${preqcdir}/twb${batch}_geno${GENO_INIT}.imiss > ${preqcdir}/twb${batch}_geno${GENO_INIT}_mind${MIND_FINAL}.indvlist

run_in_env r_env Rscript \
    01_idv_missing_plot.r $preqcdir twb${batch}_geno${GENO_INIT}.imiss

# 3. Filter to samples in the previous step, output a list of SNPs with call rate >0.98
plink2 --bfile ${datadir}/${prefix} \
    --allow-extra-chr \
    --extract ${preqcdir}/twb${batch}_geno${GENO_INIT}.snplist \
    --keep ${preqcdir}/twb${batch}_geno${GENO_INIT}_mind${MIND_FINAL}.indvlist \
    --geno $GENO_FINAL \
    --mind $MIND_FINAL \
    --make-bed \
    --write-snplist \
    --out ${preqcdir}/twb${batch}_geno${GENO_INIT}_mind${MIND_FINAL}_geno${GENO_FINAL}

# Summarize number of samples and SNPs removed at each step
qc_sum=${preqcdir}/qc_summary.txt
echo "Quality Control Summary" > $qc_sum
printf "Initial Samples: $(wc -l < ${datadir}/${prefix}.fam)\n" >> $qc_sum
printf "Initial SNPs: $(wc -l < ${datadir}/${prefix}.bim)\n" >> $qc_sum
printf "Final Samples: $(wc -l < ${preqcdir}/twb${batch}_geno${GENO_INIT}_mind${MIND_FINAL}_geno${GENO_FINAL}.fam)\n" >> $qc_sum
printf "Final SNPs: $(wc -l < ${preqcdir}/twb${batch}_geno${GENO_INIT}_mind${MIND_FINAL}_geno${GENO_FINAL}.snplist)\n" >> $qc_sum

# 4. Find out duplicated SNPs (same chr, pos, ref, alt)
plink \
    --bfile ${preqcdir}/twb${batch}_geno${GENO_INIT}_mind${MIND_FINAL}_geno${GENO_FINAL} \
    --list-duplicate-vars ids-only suppress-first \
    --out ${preqcdir}/twb${batch}_dupcheck

# 5a. Remove duplicated and monomorphic sites
# 5b. Set all the variant IDs to chr:pos:ref:alt and remove indels with length > 10
plink2 \
    --bfile ${preqcdir}/twb${batch}_geno${GENO_INIT}_mind${MIND_FINAL}_geno${GENO_FINAL} \
    --exclude ${preqcdir}/twb${batch}_dupcheck.dupvar \
    --not-chr 0 \
    --ref-from-fa --fa $REF \
    --snps-only just-acgt \
    --set-hh-missing \
    --mac 1 \
    --set-all-var-ids @:#:\$r:\$a \
    --output-chr chrM \
    --make-bed \
    --out ${preqcdir}/twb${batch}_geno${GENO_INIT}_mind${MIND_FINAL}_geno${GENO_FINAL}_dedup
