#!/usr/bin/env Rscript

# ==============================================================================
# Generate Manhattan and QQ plots from REGENIE summary statistics.
#
# Usage:
#     Rscript 04_2_regenie_plot.r <phenotype> <panel> <sv_only>
#
# Logic:
#     - Read the REGENIE output file for the given phenotype and panel.
#     - If 'sv_only' is TRUE, filter out non-SV variants.
#     - Extract additive model (ADD) results and filter P-values < 1e-5.
#     - Compute cumulative base pair positions for plotting.
#     - Generate a Manhattan plot using ggplot2.
#     - Generate a QQ plot with Genomic Control (lambda_gc) annotation.
#
# Outputs:
#     <prefix>_regenie_1e-5.txt - Table of top variants (P < 1e-5).
#     <prefix>_manhattan_plt.pdf / .jpg - Manhattan plot visualizations.
#     <prefix>_qq_plt.pdf / .jpg - QQ plot visualizations.
# ==============================================================================

args <- commandArgs(TRUE)

# Check if the correct number of arguments is provided
if (length(args) < 3) {
  stop("Usage: Rscript script.R <phenotype> <panel> <sv_only>")
}

pheno <- args[1]
panel <- args[2]  # could be "sv_snv" or "snv"
sv_only <- as.logical(args[3])

prefix <- paste(panel, "rsq0.8_maf0.005", sep = ".")  # sv_snv.rsq0.8_maf0.005
wkdir <- Sys.getenv("OUT_DIR")
if (wkdir == "") wkdir <- "/staging/biology/u4432941/sv/prs/outputs"
file <- file.path(wkdir, paste0(prefix, ".regenie_results_prs"), paste0(pheno, ".regenie"))
# /staging/biology/u4432941/sv/prs/outputs/sv_snv.rsq0.8_maf0.005.regenie_results_prs/[PHENO].regenie

library(tidyverse)

# Create output directory based on sv_only flag
if (sv_only) {
  # Filter to SV (DEL and INS) variants only
  outdir <- file.path(wkdir, paste0(prefix, ".manhattan_plot.sv_only.prs"))
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  prefix <- paste0(outdir, "/", pheno, "_")
} else {
  # Keep all variants
  outdir <- file.path(wkdir, paste0(prefix, ".manhattan_plot.all_variants.prs"))
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  prefix <- paste0(outdir, "/", pheno, "_")
}


assoc_data <- read.table(file,  header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  select(-1, -2) # Remove duplicate CHROM and GENPOS columns

# Filter assoc_data based on the bolean sv_only
if (sv_only) {
  assoc_data <- assoc_data %>%
    filter(VAR_TYPE %in% c("DEL", "INS"))
}

# Filter to include only additive model results
assoc_data <- assoc_data %>%
  filter(TEST == "ADD") %>% 
  mutate(P = 10^(-LOG10P))
assoc_data$CHROM <- factor(as.character(assoc_data$CHROM), levels = as.character(1:22))

# Output a table with P<1e-5
assoc_data %>%
  filter(P < 1e-5) %>%
  select(ID, CHROM, POS, REF, ALT, ALLELE0, ALLELE1, VAR_TYPE, N, BETA, SE, CHISQ, LOG10P, P, EXTRA) %>%
  write.table(file = paste0(prefix, "regenie_1e-5.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

# Compute cumulative positions for plotting
don <- assoc_data %>%
  group_by(CHROM) %>%
  summarise(chr_len = max(POS)) %>%
  mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
  select(-chr_len) %>%
  left_join(assoc_data, ., by = "CHROM") %>%
  arrange(CHROM, POS) %>%
  mutate(BPcum = POS + tot,
          is_highlight = ifelse(LOG10P > -log10(5e-8), "yes", "no"),
          VAR_TYPE = factor(VAR_TYPE, levels = c("SNV", "INDEL", "DEL", "INS", "Other")))

axisdf <- don %>%
  group_by(CHROM) %>%
  summarize(center = (max(as.numeric(BPcum)) + min(as.numeric(BPcum))) / 2)

# Create Manhattan plot
p <- ggplot(don, aes(x = BPcum, y = LOG10P)) +
  geom_point(aes(color = as.factor(CHROM), shape = VAR_TYPE), alpha = 0.8, size = 1.3) +
  scale_color_manual(values = rep(c("#569DC1", "#78759F"), 22)) +
  scale_shape_manual(values = c(
    SNV = 16, INDEL = 16,    # filled triangle
    DEL = 15,      # filled square
    INS = 3,       # plus
    Other = 4      # x
  )) +
  scale_x_continuous(label = axisdf$CHROM, breaks = axisdf$center) +
  scale_y_continuous(limits = c(0, NA), expand = c(0.05, 0)) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  labs(x = "Chromosome", y = "-log10(P)") +
  geom_hline(yintercept = -log10(5e-8), color = "#C84630", linetype = "dashed") +
  guides(color = "none")

ggsave(filename = paste0(prefix, "manhattan_plt.pdf"), plot = p, height = 7, width = 11.5)
ggsave(filename = paste0(prefix, "manhattan_plt.jpg"), plot = p, height = 7, width = 11.5, dpi = 500)

# QQ plot
qqplot_with_gc <- function(ps, pheno_name = "", ci = 0.95, lambda_gc = NULL) {
  ps <- na.omit(ps)
  n <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  
  gg <- ggplot(df) +
    geom_point(aes(expected, observed)) +
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_ribbon(aes(x = expected, ymin = clower, ymax = cupper), fill = "grey", alpha = 0.4) + 
    xlab(log10Pe) +
    ylab(log10Po) +
    labs(title = pheno_name)
  
  if (!is.null(lambda_gc)) {
    gg <- gg + annotate("text", x = Inf, y = -Inf, hjust = 1, vjust = -0.5,
                        label = paste("Genomic Control (lambda_gc):", round(lambda_gc, 4)))
  }
  
  return(gg)
}

# Calculate lambda_gc and add it to the plot
chisq <- qchisq(na.omit(assoc_data$P), 1, lower.tail = FALSE)
lambda_gc <- median(chisq) / qchisq(0.5, 1)

qq_plot <- qqplot_with_gc(assoc_data$P, lambda_gc = lambda_gc)

ggsave(filename = paste0(prefix, "qq_plt.pdf"), qq_plot, height = 7, width = 9)
ggsave(filename = paste0(prefix, "qq_plt.jpg"), plot = qq_plot, height = 7, width = 9, dpi = 500)

