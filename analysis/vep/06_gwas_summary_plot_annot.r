#!/usr/bin/env Rscript

# ==============================================================================
# Generate a summary plot of significant variants across all phenotypes.
#
# Usage:
#     Rscript 05_summary_plot.r [filter_to_sv] [vep_dir] [label_mode] [pleiotropy_min]
#
# Arguments:
#     filter_to_sv    : TRUE/FALSE — filter to SV variants only [default: FALSE]
#     vep_dir         : (Optional) Directory containing VEP-annotated TSV files.
#                       Files are expected to follow the naming convention:
#                       <vep_dir>/quan_sv_vep_annotated.tsv  (for Quantitative)
#                       <vep_dir>/bi_sv_vep_annotated.tsv    (for Binary)
#     label_mode      : (Optional) "significant", "pleiotropic", or "both" [default: "both"]
#     pleiotropy_min  : (Optional) Minimum pleiotropy count for pleiotropic labels [default: 5]
#
# Logic:
#     - Read all strict-filtered .regenie association results for all phenotypes.
#     - Group phenotypes into categories (e.g., Anthropometric, Cardiac, Kidney).
#     - Generate a genome-wide summary plot showing significant signals by chromosome.
#     - Optionally filter to structural variants only.
#     - If vep_dir is provided, load the type-specific VEP TSV per TYPE loop iteration
#       and overlay gene labels and vertical lines on annotated variants.
#
# Outputs:
#     <OUTPUT_DIR>/summary_plots/strict_all_sig_plot_<TYPE>.pdf
#     <OUTPUT_DIR>/summary_plots/strict_all_sig_plot_<TYPE>_annotated.pdf  (if vep_dir provided)
#     (Or strict_sv_only_plot_<TYPE>.pdf/.RData if filter_to_sv is TRUE)
# ==============================================================================

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

# Parse arguments
filter_to_sv   <- if (length(args) >= 1) as.logical(args[1]) else FALSE
vep_dir        <- if (length(args) >= 2 && args[2] != "") args[2] else NULL
label_mode     <- if (length(args) >= 3) args[3] else "both"   # "significant", "pleiotropic", or "both"
pleiotropy_min <- if (length(args) >= 4) as.integer(args[4]) else 5L

# Load annotation libraries only if vep_dir is provided
if (!is.null(vep_dir)) {
  library(ggnewscale)
}

# Debug print
cat(paste0("filter_to_sv   = ", filter_to_sv, "\n"))
cat(paste0("vep_dir        = ", ifelse(is.null(vep_dir), "(none)", vep_dir), "\n"))
cat(paste0("label_mode     = ", label_mode, "\n"))
cat(paste0("pleiotropy_min = ", pleiotropy_min, "\n"))

# Set output directory
output_dir <- Sys.getenv("OUTPUT_DIR", "/staging/biology/u4432941/sv/prs/outputs")

# Create summary_plot directory
outdir <- file.path(output_dir, "summary_plots")
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

phenotype_files_path <- file.path(output_dir, "sv_snv.rsq0.8_maf0.005.regenie_results/filtered_strict")

# Define base output path based on filter_to_sv
if (filter_to_sv) {
  plot_output_path  <- file.path(outdir, "strict_sv_only_plot.pdf")
  rdata_output_path <- file.path(outdir, "strict_sv_only_plot.RData")
} else {
  plot_output_path  <- file.path(outdir, "strict_all_sig_plot.pdf")
}

# Read the phenotype metadata to define groups and order
pheno_meta_file <- Sys.getenv("PHENO_META", "/staging/biology/u4432941/sv/prs/scripts/config/phenotype_metadata.csv")
if (!file.exists(pheno_meta_file)) {
  stop("PHENO_META file not found: ", pheno_meta_file)
}
pheno_group_df_full <- suppressMessages(read_csv(pheno_meta_file))

# ==============================================================================
# Read in phenotype association files
# ==============================================================================

all_sig_dfs <- list()

for (pheno in pheno_group_df_full$PHENOTYPE) {
  file_name <- file.path(phenotype_files_path, paste0(pheno, ".filtered.regenie"))

  if (file.exists(file_name)) {
    current_df <- suppressMessages(read_tsv(file_name, col_types = cols(.default = "c"))) %>%
      select(-1, -2) %>%
      mutate(PHENOTYPE = pheno)

    if (filter_to_sv) {
      current_df <- current_df %>%
        filter(VAR_TYPE %in% c("DEL", "INS"))
    }

    cat(sprintf("Phenotype %s: %d significant variants retained for plotting\n", pheno, nrow(current_df)))
    all_sig_dfs[[pheno]] <- current_df
  } else {
    cat(sprintf("Phenotype %s: No significant results file found (0 variants)\n", pheno))
  }
}

# Combine all significant results and infer column types
sig_df <- suppressMessages(bind_rows(all_sig_dfs) %>% type_convert())

# ==============================================================================
# Loop over each phenotype TYPE (Binary / Quantitative)
# ==============================================================================

types <- unique(pheno_group_df_full$TYPE)

for (current_type in types) {
  cat(paste0("\n=== Generating plot for TYPE: ", current_type, " ===\n"))

  # --------------------------------------------------------------------------
  # Load type-specific VEP annotation TSV
  # --------------------------------------------------------------------------
  vep_df <- NULL

  if (!is.null(vep_dir)) {
    # Map TYPE to file name prefix: Quantitative -> quan, Binary -> bi
    type_prefix <- switch(current_type,
      "Quantitative" = "quan",
      "Binary"       = "bi",
      tolower(current_type)  # fallback: other TYPE values are lowercased directly
    )
    vep_tsv <- file.path(vep_dir, paste0(type_prefix, "_sv_vep_annotated.tsv"))

    if (file.exists(vep_tsv)) {
      message("Loading VEP annotation from: ", vep_tsv)
      vep_df <- read.delim(vep_tsv, stringsAsFactors = FALSE) %>%
        mutate(
          pleiotropy_count   = str_count(associated_phenotypes, ",") + 1,
          is_top_pleiotropic = pleiotropy_count >= pleiotropy_min
        ) %>%
        select(all_of(c("ID", "associated_phenotypes", "Consequence", "IMPACT",
                        "SYMBOL", "pleiotropy_count", "is_top_pleiotropic")))
    } else {
      message("VEP TSV not found for TYPE '", current_type, "': ", vep_tsv, " — skipping annotation.")
    }
  }

  # --------------------------------------------------------------------------
  # Filter metadata and results for this TYPE
  # --------------------------------------------------------------------------
  pheno_group_df <- pheno_group_df_full %>%
    filter(TYPE == current_type) %>%
    mutate(PHENO_INDEX = rev(seq_along(PHENOTYPE)))

  type_sig_df <- sig_df %>% filter(PHENOTYPE %in% pheno_group_df$PHENOTYPE)

  if (nrow(type_sig_df) == 0) {
    cat(paste0("No significant results for TYPE: ", current_type, ". Skipping...\n"))
    next
  }

  df <- left_join(type_sig_df, pheno_group_df, by = "PHENOTYPE")

  df$CHROM     <- factor(df$CHROM,     levels = as.character(c(1:22)))
  df$GROUP     <- factor(df$GROUP,     levels = unique(pheno_group_df$GROUP))
  df$PHENOTYPE <- factor(df$PHENOTYPE, levels = pheno_group_df$PHENOTYPE)
  df$VAR_TYPE  <- factor(df$VAR_TYPE,  levels = c("SNV", "INDEL", "DEL", "INS"))

  num_groups <- length(unique(pheno_group_df$GROUP))
  my_colors  <- colorRampPalette(RColorBrewer::brewer.pal(9, "Pastel1"))(num_groups)

  group_bands <- pheno_group_df %>%
    group_by(GROUP) %>%
    summarize(ymin = min(PHENO_INDEX) - 0.5,
              ymax = max(PHENO_INDEX) + 0.5) %>%
    mutate(xmin = -Inf, xmax = Inf)

  # --------------------------------------------------------------------------
  # Base summary plot
  # --------------------------------------------------------------------------
  plt <- ggplot(df, aes(x = POS, y = PHENO_INDEX)) +
    geom_rect(data = group_bands,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = GROUP),
              alpha = 0.2, color = NA, inherit.aes = FALSE) +
    geom_point(aes(shape = VAR_TYPE), size = 2, stroke = 0.6) +
    scale_x_continuous(
      labels = NULL,
      name   = NULL,
      expand = expansion(mult = c(0.001, 0.001))
    ) +
    scale_y_continuous(
      breaks = pheno_group_df$PHENO_INDEX,
      labels = pheno_group_df$PHENOTYPE,
      name   = "Phenotype"
    ) +
    scale_shape_manual(values = c(SNV   = 16, INDEL = 16,
                                  DEL   = 15,
                                  INS   = 3,
                                  Other = 4)) +
    scale_fill_manual(values = my_colors) +
    facet_grid(
      . ~ CHROM,
      scales = "free_x",
      space  = "free_x",
      switch = "x"
    ) +
    coord_cartesian(clip = "off") +
    theme_minimal() +
    theme(
      panel.grid       = element_blank(),
      panel.spacing.x  = unit(0.1, "lines"),
      strip.background = element_rect(fill = "gray95", color = "black"),
      strip.text       = element_text(size = 9, face = "bold"),
      axis.text.y      = element_text(size = 8),
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      legend.position  = "right",
      axis.title       = element_blank(),
      plot.margin      = margin(10, 10, 10, 10)
    ) +
    labs(
      shape = "Variant Type",
      fill  = "Phenotype Group"
    )

  # Define type-specific output paths
  type_plot_output_path <- sub("\\.pdf$", paste0("_", current_type, ".pdf"), plot_output_path)

  if (filter_to_sv) {
    type_rdata_output_path <- sub("\\.RData$", paste0("_", current_type, ".RData"), rdata_output_path)
    save(plt, df, group_bands, file = type_rdata_output_path)
    cat(paste0("R plot object saved to: ", type_rdata_output_path, "\n"))
  } else {
    cat("No R plot object saved as filter_to_sv is FALSE.\n")
  }

  num_phenos  <- nrow(pheno_group_df)
  plot_height <- max(6, num_phenos * 0.25)

  pdf(type_plot_output_path, width = 15, height = plot_height)
  print(plt)
  dev.off()
  cat(paste0("Plot saved to: ", type_plot_output_path, "\n"))

  # --------------------------------------------------------------------------
  # Annotated plot — only if VEP TSV was successfully loaded for this TYPE
  # --------------------------------------------------------------------------
  if (!is.null(vep_df)) {

    df_annot <- df %>%
      left_join(vep_df, by = "ID")

    max_pheno_index <- max(pheno_group_df$PHENO_INDEX)

    # Label set 1: all significant variants with a gene symbol
    gene_label_significant <- df_annot %>%
      filter(!is.na(SYMBOL), SYMBOL != ".") %>%
      distinct(SYMBOL, CHROM, .keep_all = TRUE)

    # Label set 2: pleiotropic variants with a gene symbol
    gene_label_pleiotropic <- df_annot %>%
      filter(is_top_pleiotropic == TRUE, !is.na(SYMBOL), SYMBOL != ".") %>%
      distinct(SYMBOL, CHROM, .keep_all = TRUE)

    label_df <- switch(label_mode,
      "significant" = gene_label_significant,
      "pleiotropic" = gene_label_pleiotropic,
      "both"        = bind_rows(gene_label_significant, gene_label_pleiotropic) %>%
                        distinct(SYMBOL, CHROM, .keep_all = TRUE),
      stop("label_mode must be one of: 'significant', 'pleiotropic', 'both'")
    )

    annot_plt <- plt

    if (nrow(label_df) > 0) {
      annot_plt <- annot_plt +
        # Vertical guide lines at each labelled variant (correct within each CHROM facet)
        geom_vline(data = label_df, aes(xintercept = POS),
                   linewidth = 0.1, color = "grey40", linetype = "dotted",
                   inherit.aes = FALSE) +

        # Reset fill scale so gene labels use their own fill independent of group bands
        new_scale_fill() +

        # Gene name labels above the top phenotype row, colour-coded by IMPACT
        geom_label(
          data        = label_df,
          aes(x     = POS,
              y     = max_pheno_index + 1.5,
              label = SYMBOL,
              fill  = ifelse(IMPACT == "HIGH", "high", "normal")),
          fontface    = "italic",
          color       = "black",
          label.size  = 0.2,
          angle       = 90,
          hjust       = 0,
          size        = 3,
          inherit.aes = FALSE,
          show.legend = FALSE
        ) +
        scale_fill_manual(values = c("high" = "#fcbcb2", "normal" = "white"))
    } else {
      message("No variants passed the annotation filter for label_mode = '", label_mode,
              "' (TYPE: ", current_type, "). Annotated plot will be identical to the base plot.")
    }

    type_annot_output_path <- sub("\\.pdf$", paste0("_", current_type, "_annotated.pdf"), plot_output_path)

    # Extra height to accommodate rotated gene labels
    pdf(type_annot_output_path, width = 15, height = plot_height + 1.5)
    print(annot_plt)
    dev.off()
    cat(paste0("Annotated plot saved to: ", type_annot_output_path, "\n"))
  }
}
