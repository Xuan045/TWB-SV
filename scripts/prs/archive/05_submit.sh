#!/usr/bin/bash

# merge pgen
SV_MERGE=$(sbatch --parsable --export=ALL,TARGET=sv_snv 05_merge_pgen.sh)
SNV_MERGE=$(sbatch --parsable --export=ALL,TARGET=snv 05_merge_pgen.sh)

# sv_snv panel
SV_QUAN=$(sbatch --parsable --dependency=afterany:$SV_MERGE --export=ALL,TARGET=sv_snv 05_2_run_prsice_quan.sh)
SV_BI=$(sbatch --parsable --dependency=afterany:$SV_MERGE --export=ALL,TARGET=sv_snv 05_3_run_prsice_bi.sh)

# snv panel
SNV_QUAN=$(sbatch --parsable --dependency=afterany:$SNV_MERGE --export=ALL,TARGET=snv 05_2_run_prsice_quan.sh)
SNV_BI=$(sbatch --parsable --dependency=afterany:$SNV_MERGE --export=ALL,TARGET=snv 05_3_run_prsice_bi.sh)
