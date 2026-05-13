# ==============================================================================
# PRS performance visualization
#
# Part 1: Per-phenotype violin + boxplot with paired lines
# Part 2: Combined bar plot with 95% CI across all phenotypes
#
# Input:
#     <data_dir>/<PHENO>_cv_results.tsv  : per-repeat R² for both panels
# Output:
#     <out_dir>/<PHENO>_violin.pdf       : per-phenotype plot
#     <out_dir>/combined_barplot.pdf     : all phenotypes combined
# ==============================================================================

library(tidyverse)
library(ggpubr)

wkdir <- Sys.getenv("OUT_DIR")
if (wkdir == "") wkdir <- "/staging/biology/u4432941/sv/prs/outputs"

data_dir <- file.path(wkdir, "prsice", "sv_snv.quan")
out_dir  <- file.path(wkdir, "prsice", "comparison_plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Phenotype list
pheno_list <- c(
  "BODY_HEIGHT", "BODY_WEIGHT", "BMI", "BODY_FAT_RATE",
  "BODY_WAISTLINE", "BODY_BUTTOCKS", "WHR",
  "T_SCORE", "Z_SCORE",
  "SYSTOLIC_PRESSURE", "DIASTOLIC_PRESSURE", "HEARTBEAT_SPEED",
  "RBC", "WBC", "HB", "HCT", "PLATELET",
  "HBA1C", "AC_GLUCOSE",
  "T_CHO", "TG", "HDL_C", "LDL_C",
  "T_BILIRUBIN", "ALBUMIN", "SGOT", "SGPT", "GAMMA_GT", "AFP",
  "BUN", "CREATININE", "URIC_ACID", "EGFRcr"
)

# Colours: SNV+SV = coral, SNV = steel blue
panel_colors <- c("SNV+SV" = "#E07B6A", "SNV" = "#5B8DB8")

# ==============================================================================
# Helper: prepare long-format data for one phenotype
# ==============================================================================

load_plot_df <- function(pheno, data_dir) {
  f <- file.path(data_dir, paste0(pheno, "_cv_results.tsv"))
  if (!file.exists(f)) {
    warning("File not found, skipping: ", f)
    return(NULL)
  }
  read.table(f, header = TRUE) %>%
    select(repeat_id, t1_test_r2, t2_test_r2) %>%
    pivot_longer(
      cols      = c(t1_test_r2, t2_test_r2),
      names_to  = "model",
      values_to = "r2"
    ) %>%
    mutate(
      model     = recode(model, t1_test_r2 = "SNV+SV", t2_test_r2 = "SNV"),
      phenotype = pheno
    )
}

# ==============================================================================
# Script 1: Per-phenotype violin + boxplot
# ==============================================================================

plot_single_pheno <- function(pheno, data_dir, out_dir) {
  df <- read.table(file.path(data_dir, paste0(pheno, "_cv_results.tsv")), header = TRUE)
  plot_df <- load_plot_df(pheno, data_dir)
  if (is.null(plot_df)) return(invisible(NULL))

  # Paired t-test: SNV+SV vs SNV
  stat.test <- t.test(df$t1_test_r2, df$t2_test_r2, paired = TRUE)
  p_label   <- paste0("paired p = ", signif(stat.test$p.value, 3))

  p <- ggplot(plot_df, aes(x = model, y = r2, fill = model, color = model)) +
    geom_violin(trim = FALSE, alpha = 0.25, linewidth = 0.4) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7, linewidth = 0.5, color = "grey30") +
    # Paired lines
    geom_line(aes(group = repeat_id), alpha = 0.12, color = "grey50", linewidth = 0.3) +
    geom_point(alpha = 0.35, size = 1.2) +
    # p-value label
    annotate("text", x = 1.5, y = max(plot_df$r2) * 1.02,
             label = p_label, size = 3.5, color = "grey30") +
    scale_fill_manual(values  = panel_colors) +
    scale_color_manual(values = panel_colors) +
    scale_x_discrete(limits = c("SNV+SV", "SNV")) +
    labs(
      x     = NULL,
      y     = expression(Incremental ~ R^2),
      title = pheno
    ) +
    theme_bw(base_size = 13) +
    theme(
      legend.position  = "none",
      panel.grid.major.x = element_blank(),
      plot.title       = element_text(face = "bold", hjust = 0.5)
    )

  ggsave(file.path(out_dir, paste0(pheno, "_violin.pdf")),
         plot = p, width = 4, height = 5)
  ggsave(file.path(out_dir, paste0(pheno, "_violin.jpg")),
         plot = p, width = 4, height = 5, dpi = 300)

  invisible(p)
}

# Run for all phenotypes
walk(pheno_list, plot_single_pheno, data_dir = data_dir, out_dir = out_dir)

# ==============================================================================
# Combined horizontal bar plot with 95% CI + delta R² annotation
# ==============================================================================

# Load and combine all phenotypes — missing from current script
all_df <- map(pheno_list, load_plot_df, data_dir = data_dir) %>%
  compact() %>%
  bind_rows() %>%
  # Make SNV+SV appear before SNV in legend and plots
  mutate(model = factor(model, levels = c("SNV+SV", "SNV")))

# Compute mean and 95% CI per phenotype x model
summary_df <- all_df %>%
  group_by(phenotype, model) %>%
  summarise(
    mean_r2  = mean(r2),
    ci_lower = quantile(r2, 0.025),
    ci_upper = quantile(r2, 0.975),
    .groups  = "drop"
  )

# Compute paired t-test p-value + delta R² per phenotype
pval_df <- map_dfr(pheno_list, function(pheno) {
  f <- file.path(data_dir, paste0(pheno, "_cv_results.tsv"))
  if (!file.exists(f)) return(NULL)
  df <- read.table(f, header = TRUE)
  tt    <- t.test(df$t1_test_r2, df$t2_test_r2, paired = TRUE)
  delta <- mean(df$t1_test_r2 - df$t2_test_r2)  # mean (SNV+SV) - (SNV)
  tibble(
    phenotype = pheno,
    p_value   = tt$p.value,
    delta_r2  = delta,
    # x position: just beyond the rightmost CI for this phenotype
    x_pos = max(
      summary_df$ci_upper[summary_df$phenotype == pheno],
      na.rm = TRUE
    )
  )
}) %>%
  mutate(
    p_label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    ),
    # Delta label: show sign explicitly, e.g. "+0.0023**" or "-0.0001 ns"
    delta_label = paste0(
      ifelse(delta_r2 >= 0, "+", ""),
      formatC(delta_r2, format = "f", digits = 4),
      " ", p_label
    ),
    delta_color = ifelse(delta_r2 >= 0, "#E07B6A", "#5B8DB8")
  )

# Order phenotypes by SNV+SV mean R² (ascending for horizontal plot — top = highest)
pheno_order <- summary_df %>%
  filter(model == "SNV+SV") %>%
  arrange(mean_r2) %>%
  pull(phenotype)

summary_df <- summary_df %>%
  mutate(phenotype = factor(phenotype, levels = pheno_order))
pval_df <- pval_df %>%
  mutate(phenotype = factor(phenotype, levels = pheno_order))

# x-axis upper limit: leave room for delta label on the right
x_max <- max(summary_df$ci_upper, na.rm = TRUE) * 1.35

combined_plt <- ggplot(summary_df,
                       aes(x = mean_r2, y = phenotype,
                           fill = model, color = model)) +
  geom_col(
    position = position_dodge(width = 0.7),
    width    = 0.6,
    alpha    = 0.85
  ) +
  geom_errorbar(
    aes(xmin = ci_lower, xmax = ci_upper),
    position  = position_dodge(width = 0.7),
    width     = 0.25,
    linewidth = 0.5,
    color     = "grey30"
  ) +
  # Delta R² label: (SNV+SV) - (SNV-only), with significance stars
  geom_text(
    data        = pval_df,
    aes(x = x_pos * 1.05, y = phenotype, label = delta_label),
    inherit.aes = FALSE,
    hjust       = 0,
    size        = 3.2,
    color       = pval_df$delta_color,
    fontface    = "plain"
  ) +
  # Vertical reference line at x = 0
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey50") +
  scale_fill_manual(values  = panel_colors, name = "Panel") +
  scale_color_manual(values = panel_colors, name = "Panel") +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0)),
    limits = c(0, x_max)
  ) +
  labs(
    y     = NULL,
    x     = expression(Incremental ~ R^2 ~ "(mean ± 95% CI)"),
    title = "PRS performance: SNV+SV vs SNV-only"
  ) +
  # Add a secondary axis label explaining the delta annotation
  annotate("text",
           x = x_max, y = 0.3,
           label = expression(Delta * R^2 ~ "(SNV+SV \u2212 SNV-only)"),
           hjust = 1, vjust = 0, size = 3, color = "grey40",
           fontface = "italic") +
  theme_bw(base_size = 13) +
  theme(
    axis.text.y      = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    legend.position  = "top",
    plot.title       = element_text(face = "bold", hjust = 0.5)
  )

# Height scales with number of phenotypes
plot_height <- max(5, length(pheno_list) * 0.35)

ggsave(file.path(out_dir, "combined_barplot.pdf"),
       plot   = combined_plt,
       width  = 9,
       height = plot_height)
ggsave(file.path(out_dir, "combined_barplot.jpg"),
       plot   = combined_plt,
       width  = 9,
       height = plot_height,
       dpi    = 300)

cat("Done. Horizontal plot saved to:", out_dir, "\n")
