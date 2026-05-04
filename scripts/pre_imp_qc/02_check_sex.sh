#!/usr/bin/bash
#SBATCH -A MST109178        # Account name/project number
#SBATCH -J check_sex         # Job name
#SBATCH -p ngs92G           # Partition Name 等同PBS裡面的 -q Queue name
#SBATCH -c 14               # 使用的core數 請參考Queue資源設定
#SBATCH --mem=92g           # 使用的記憶體量 請參考Queue資源設定
#SBATCH -o /dev/null
#SBATCH -e /dev/null 

# ---- Source Configuration ----
source "$(dirname "$0")/../config.sh"

# ---- Environment ----
# Environment is initialized in config.sh
activate_env assoc_env

# ---- Parameters ----
GENO_INIT=0.05
GENO_FINAL=0.02
MIND_FINAL=0.02
batch='2'

# ---- Redirect LOG file ----
logfile=${LOG_DIR}/02_check_sex_v${batch}.log
exec > "$logfile" 2>&1

preqcdir=${OUT_DIR}/pre_imp_qc_v${batch}
BFILE=${preqcdir}/twb${batch}_geno${GENO_INIT}_mind${MIND_FINAL}_geno${GENO_FINAL}_dedup
PARA=${preqcdir}/twb${batch}

# Perform LD pruning: chrX
plink \
    --bfile ${BFILE} \
    --chr X,Y,XY \
    --maf 0.05 \
    --snps-only just-acgt \
    --indep-pairwise 200 100 0.1 \
    --out ${PARA}.chrX

# Check sex
plink \
    --bfile ${BFILE} \
    --extract ${PARA}.chrX.prune.in \
    --check-sex \
    --out ${PARA}.chrX

# Plot
run_in_env r_env \
    Rscript 02_check_sex.r \
    $preqcdir $PARA 0.2
