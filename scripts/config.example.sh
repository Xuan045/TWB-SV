#!/usr/bin/bash
# -----------------------------------------------------------------------------
# Central Configuration for PRS Pipeline
# -----------------------------------------------------------------------------

# Base Project Directory
# TODO: Update this to point to the root directory where you cloned the repository
export PROJECT_DIR="/path/to/your/project/dir"

# -----------------------------------------------------------------------------
# Environment Management Setup
# -----------------------------------------------------------------------------
# Supported options: "micromamba", "conda"
export ENV_MANAGER="micromamba"

# Initialize base environment based on the selected manager
if [ "$ENV_MANAGER" == "micromamba" ]; then
    # TODO: Replace with your micromamba path
    export MAMBA_EXE="/path/to/bin/micromamba"
    export MAMBA_ROOT_PREFIX="/path/to/micromamba/root"
    eval "$("$MAMBA_EXE" shell hook --shell bash --prefix "$MAMBA_ROOT_PREFIX")"
elif [ "$ENV_MANAGER" == "conda" ]; then
    # TODO: Load your conda module here if on an HPC, or comment out if not needed
    ml biology
    ml Anaconda/Anaconda3
fi

# Function to activate an environment
activate_env() {
    local env_name=$1
    if [ "$ENV_MANAGER" == "micromamba" ]; then
        micromamba activate "$env_name"
    elif [ "$ENV_MANAGER" == "conda" ]; then
        conda activate "$env_name" 
    fi
}

# Function to run a command inside an environment
run_in_env() {
    local env_name=$1
    shift # Remove env_name from arguments
    if [ "$ENV_MANAGER" == "micromamba" ]; then
        micromamba run -n "$env_name" "$@"
    elif [ "$ENV_MANAGER" == "conda" ]; then
        conda run -n "$env_name" "$@"
    fi
}

# -----------------------------------------------------------------------------
# Core Directories
# -----------------------------------------------------------------------------

export OUT_DIR="${PROJECT_DIR}/outputs"
export LOG_DIR="${PROJECT_DIR}/logs/prs"
export PRSICE_LOG_DIR="${PROJECT_DIR}/logs/prsice"

# External / Shared Resource Paths
# TODO: Provide the correct paths for your reference data
export LAB_INFO="/path/to/your/lab_info.csv"
export SURVEY_DATA="/path/to/your/full_df_survey.csv"
export LAB_DATA="/path/to/your/full_df_lab.csv"
export PHENO_FILE="${PROJECT_DIR}/pheno_summary/twb2_pheno_covar.txt"
export PHENO_META="${PROJECT_DIR}/scripts/config/phenotype_metadata.csv"

# Function to extract comma-separated phenotype lists by type (Quantitative or Binary)
get_pheno_list() {
    local ptype=$1
    awk -F, -v type="$ptype" 'NR>1 && $2==type {printf "%s,", $1}' "$PHENO_META" | sed 's/,$//'
}

# Function to safely exclude a phenotype from a comma-separated list
exclude_pheno() {
    local list="$1"
    local to_exclude="$2"
    echo "$list" | tr ',' '\n' | grep -vFx "$to_exclude" | paste -sd, -
}

# External Tool Paths
export PRSICE_BIN="/path/to/your/PRSice_linux"
export SHAPEIT5_DIR="/path/to/your/shapeit5_v5.1.1"
export MINIMAC4_BIN="/path/to/your/minimac4"

# Pre-computed Data Paths
export BCF_PHASED_V2="${OUT_DIR}/imputation_v2/phasing/phased.chr1.bcf"
export BFILE_PCA="${OUT_DIR}/pca_v2/twb2_preqc_sexcheck_hetcheck_EAS"
export KIN0_FILE="${OUT_DIR}/check_ibd/twb2.ldpr.2degree.kin0"

# Reference Files
export FASTA_REF="/path/to/your/GRCh38_full_analysis_set_plus_decoy_hla.fa"
export CHAIN_FILE="${PROJECT_DIR}/reference/hg19ToHg38.over.chain"
export LD_EXCLUDE="${PROJECT_DIR}/reference/long_range_LD_intervals.txt"
export GMAP_PREFIX="/path/to/your/gmap.b38/chr"
export IMP_PANEL_PREFIX="/path/to/your/imp_panel"
