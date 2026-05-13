#!/usr/bin/env Rscript

# ==============================================================================
# 100-repeat validation/testing split for C+T PRS evaluation.
# Compares sv_snv vs snv by incremental R² against a covariate-only null model.
#
# Usage:
#     Rscript 05_2_prsice_cv.r <phenotype> <target1> <target2> [n_repeat] [n_validation]
#
# Arguments:
#     phenotype    : Phenotype name
#     target1      : First panel  (e.g., sv_snv)
#     target2      : Second panel (e.g., snv)
#     n_repeat     : Number of repeats [default: 100]
#     n_validation : Validation set size [default: 5000]
#
# Logic:
#     For each repeat:
#       1. Use the same random validation/testing split for both targets.
#       2. In validation: select the threshold with highest incremental R²
#          (vs covariate-only null) independently for each target.
#       3. In testing: evaluate incremental R² for each target using its
#          own selected threshold.
#     Output per-repeat results and a summary table for both targets.
#
# Outputs:
#     <outdir>/<phenotype>_cv_results.tsv   : per-repeat R² for both targets
#     <outdir>/<phenotype>_cv_summary.tsv   : mean / SD / 95% CI per target
# ==============================================================================

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript 05_2_prsice_cv.r <phenotype> <target1> <target2> [n_repeat] [n_validation]")
}

pheno        <- args[1]
target1      <- args[2]   # sv_snv
target2      <- args[3]   # snv
n_repeat     <- if (length(args) >= 4) as.integer(args[4]) else 100L
n_validation <- if (length(args) >= 5) as.integer(args[5]) else 5000L

output_dir <- Sys.getenv("OUT_DIR", "/staging/biology/u4432941/sv/prs/outputs")
pheno_file <- Sys.getenv("PHENO_FILE")

cat(sprintf("Phenotype    : %s\n", pheno))
cat(sprintf("Target 1     : %s\n", target1))
cat(sprintf("Target 2     : %s\n", target2))
cat(sprintf("N repeats    : %d\n", n_repeat))
cat(sprintf("N validation : %d\n", n_validation))

# ==============================================================================
# Helper functions
# ==============================================================================

load_score <- function(target, pheno, output_dir) {
  score_file <- file.path(output_dir, "prsice", paste0(target, ".quan"),
                          paste0(pheno, ".all_score"))
  if (!file.exists(score_file)) stop("Score file not found: ", score_file)
  df <- read.table(score_file, header = TRUE, stringsAsFactors = FALSE)
  colnames(df)[1:2] <- c("FID", "IID")
  df
}

incremental_r2 <- function(df, prs_col, pheno, cov_formula) {
  full_formula <- update(cov_formula, paste(". ~ . +", prs_col))
  null_r2 <- summary(lm(cov_formula,  data = df))$r.squared
  full_r2 <- summary(lm(full_formula, data = df))$r.squared
  full_r2 - null_r2
}

select_threshold_and_eval <- function(val_df, test_df, threshold_cols, pheno, cov_formula) {
  # Select best threshold on validation set
  val_r2 <- sapply(threshold_cols, function(col) {
    incremental_r2(val_df, col, pheno, cov_formula)
  })
  best_col    <- threshold_cols[which.max(val_r2)]
  best_val_r2 <- max(val_r2)

  # Evaluate on testing set
  test_r2 <- incremental_r2(test_df, best_col, pheno, cov_formula)

  list(best_threshold = best_col, val_r2 = best_val_r2, test_r2 = test_r2)
}

# ==============================================================================
# Load data
# ==============================================================================

if (!file.exists(pheno_file)) stop("PHENO_FILE not found: ", pheno_file)
pheno_df <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)

# Load scores for both targets
score1 <- load_score(target1, pheno, output_dir)
score2 <- load_score(target2, pheno, output_dir)

threshold_cols1 <- setdiff(colnames(score1), c("FID", "IID"))
threshold_cols2 <- setdiff(colnames(score2), c("FID", "IID"))

# Merge each target's scores with phenotype/covariates
# Suffix scores to avoid column name collisions
data1 <- inner_join(score1, pheno_df, by = "IID") %>%
  filter(!is.na(.data[[pheno]]))
data2 <- inner_join(score2, pheno_df, by = "IID") %>%
  filter(!is.na(.data[[pheno]]))

# Use the intersection of IIDs so both targets use the same individuals
common_ids <- intersect(data1$IID, data2$IID)
data1 <- data1 %>% filter(IID %in% common_ids)
data2 <- data2 %>% filter(IID %in% common_ids)

cat(sprintf("Individuals used (intersection): %d\n", length(common_ids)))
cat(sprintf("Thresholds | %s: %d | %s: %d\n",
            target1, length(threshold_cols1),
            target2, length(threshold_cols2)))

# Covariate-only null model formula
cov_formula <- as.formula(paste(pheno, "~ AGE + SEX +",
                                paste(sprintf("PC%d", 1:10), collapse = " + ")))

# ==============================================================================
# 100-repeat validation / testing split
# ==============================================================================

set.seed(42)
results <- vector("list", n_repeat)

for (i in seq_len(n_repeat)) {
  if (i %% 10 == 0) cat(sprintf("  Repeat %d / %d\n", i, n_repeat))

  # 1. Same random split applied to both targets
  val_ids  <- sample(common_ids, size = n_validation, replace = FALSE)
  test_ids <- setdiff(common_ids, val_ids)

  val1  <- data1 %>% filter(IID %in% val_ids)
  test1 <- data1 %>% filter(IID %in% test_ids)
  val2  <- data2 %>% filter(IID %in% val_ids)
  test2 <- data2 %>% filter(IID %in% test_ids)

  # 2 & 3. Select threshold on validation, evaluate on testing — independently per target
  res1 <- select_threshold_and_eval(val1, test1, threshold_cols1, pheno, cov_formula)
  res2 <- select_threshold_and_eval(val2, test2, threshold_cols2, pheno, cov_formula)

  results[[i]] <- tibble(
    repeat_id          = i,
    n_val              = length(val_ids),
    n_test             = length(test_ids),
    # target1 (sv_snv)
    t1_best_threshold  = res1$best_threshold,
    t1_val_r2          = res1$val_r2,
    t1_test_r2         = res1$test_r2,
    # target2 (snv)
    t2_best_threshold  = res2$best_threshold,
    t2_val_r2          = res2$val_r2,
    t2_test_r2         = res2$test_r2
  )
}

results_df <- bind_rows(results)

# ==============================================================================
# Summary: mean / SD / 95% CI per target
# ==============================================================================

summarise_r2 <- function(x) {
  tibble(
    mean_r2  = mean(x),
    sd_r2    = sd(x),
    ci_lower = quantile(x, 0.025),
    ci_upper = quantile(x, 0.975)
  )
}

summary_df <- bind_rows(
  summarise_r2(results_df$t1_test_r2) %>% mutate(target = target1, phenotype = pheno),
  summarise_r2(results_df$t2_test_r2) %>% mutate(target = target2, phenotype = pheno)
) %>%
  select(phenotype, target, mean_r2, sd_r2, ci_lower, ci_upper)

# ==============================================================================
# Save outputs
# ==============================================================================

# Use target1's output dir as the shared output location
cv_outdir <- file.path(output_dir, "prsice", paste0(target1, ".quan"))

cv_results_path <- file.path(cv_outdir, paste0(pheno, "_cv_results.tsv"))
cv_summary_path <- file.path(cv_outdir, paste0(pheno, "_cv_summary.tsv"))

write.table(results_df, file = cv_results_path, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(summary_df, file = cv_summary_path, sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("\nDone. Results saved to:\n  %s\n  %s\n", cv_results_path, cv_summary_path))
cat(sprintf("\n%-10s | mean incr.R² | SD     | 95%% CI\n", "Target"))
cat(sprintf("%-10s | %.4f       | %.4f | [%.4f, %.4f]\n",
            target1,
            summary_df$mean_r2[summary_df$target == target1],
            summary_df$sd_r2[summary_df$target == target1],
            summary_df$ci_lower[summary_df$target == target1],
            summary_df$ci_upper[summary_df$target == target1]))
cat(sprintf("%-10s | %.4f       | %.4f | [%.4f, %.4f]\n",
            target2,
            summary_df$mean_r2[summary_df$target == target2],
            summary_df$sd_r2[summary_df$target == target2],
            summary_df$ci_lower[summary_df$target == target2],
            summary_df$ci_upper[summary_df$target == target2]))
