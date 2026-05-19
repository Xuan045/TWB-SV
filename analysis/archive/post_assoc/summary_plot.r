# ==============================================================================
# Generate a summary plot of significant variants across all phenotypes.
#
# Usage:
#     Rscript 05_summary_plot.r [filter_to_sv]
#
# Logic:
#     - Read all strict-filtered .regenie association results for all phenotypes.
#     - Group phenotypes into categories (e.g., Anthropometric, Cardiac, Kidney).
#     - Generate a genome-wide summary plot showing significant signals by chromosome.
#     - Optionally filter to structural variants only.
#
# Outputs:
#     <OUTPUT_DIR>/sv_snv.rsq0.8_maf0.005.regenie_results/summary_plots/strict_all_sig_plot.pdf
#     (Or strict_sv_only_plot.pdf if filter_to_sv is TRUE)
# ==============================================================================

library(tidyverse)

# Set output directory
output_dir <- Sys.getenv("OUTPUT_DIR", "/staging/biology/u4432941/sv/prs/outputs")

# Create summary_plot directory
outdir <- file.path(output_dir, "summary_plots")
if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
}

phenotype_files_path <- file.path(output_dir, "sv_snv.rsq0.8_maf0.005.regenie_results/filtered_strict")

args <- commandArgs(trailingOnly = TRUE)

# If filter to SV
# Default to FALSE if not provided
filter_to_sv <- if (length(args) >= 1) as.logical(args[1]) else FALSE

# Debug print
cat(paste0("filter_to_sv = ", filter_to_sv, "\n"))

# Define output path based on argument
if (filter_to_sv) {
  plot_output_path <- file.path(outdir, "strict_sv_only_plot.pdf")
  rdata_output_path <- file.path(outdir, "strict_sv_only_plot.RData")
} else {
  plot_output_path <- file.path(outdir, "strict_all_sig_plot.pdf")
}

# Read the phenotype metadata to define groups and order
pheno_meta_file <- Sys.getenv("PHENO_META", "/staging/biology/u4432941/sv/prs/scripts/config/phenotype_metadata.csv")
if (!file.exists(pheno_meta_file)) {
  stop("PHENO_META file not found: ", pheno_meta_file)
}
pheno_group_df_full <- suppressMessages(read_csv(pheno_meta_file))

# --- Read in phenotype files ---
# Initialize an empty list to store significant results from each phenotype
all_sig_dfs <- list()

# Loop through each phenotype to read its file and filter for significant variants
for (pheno in pheno_group_df_full$PHENOTYPE) {
  file_name <- file.path(phenotype_files_path, paste0(pheno, ".filtered.regenie"))

  if (file.exists(file_name)) {
    # Read the file using read_tsv for tab-separated files (common for GWAS results)
    # Ensure column names match your file exactly, especially 'LOG10P', 'CHROM', 'VAR_TYPE'
    current_df <- suppressMessages(read_tsv(file_name, col_types = cols(.default = "c"))) %>%
      select(-1, -2) %>%
      mutate(PHENOTYPE = pheno)

    # Filter for SV if specified
    if (filter_to_sv) {
      current_df <- current_df %>%
        filter(VAR_TYPE %in% c("DEL", "INS"))
    }
    
    cat(sprintf("Phenotype %s: %d significant variants retained for plotting\n", pheno, nrow(current_df)))

    # Store the significant results
    all_sig_dfs[[pheno]] <- current_df
  } else {
    # No warning if file doesn't exist, just a log message
    cat(sprintf("Phenotype %s: No significant results file found (0 variants)\n", pheno))
  }
}

# Combine all significant results into a single data frame and guess types
sig_df <- suppressMessages(bind_rows(all_sig_dfs) %>% type_convert())


# Loop over each phenotype TYPE (Binary / Quantitative)
types <- unique(pheno_group_df_full$TYPE)

for (current_type in types) {
  cat(paste0("\n=== Generating plot for TYPE: ", current_type, " ===\n"))
  
  # Filter metadata for this type
  pheno_group_df <- pheno_group_df_full %>% 
    filter(TYPE == current_type) %>%
    mutate(PHENO_INDEX = rev(seq_along(PHENOTYPE)))
  
  # Filter significant results for this type
  type_sig_df <- sig_df %>% filter(PHENOTYPE %in% pheno_group_df$PHENOTYPE)
  
  if (nrow(type_sig_df) == 0) {
    cat(paste0("No significant results for TYPE: ", current_type, ". Skipping...\n"))
    next
  }

  df <- left_join(type_sig_df, pheno_group_df, by = "PHENOTYPE")
  
  df$CHROM <- factor(df$CHROM, levels = as.character(c(1:22)))
  df$GROUP <- factor(df$GROUP, levels = unique(pheno_group_df$GROUP))
  df$PHENOTYPE <- factor(df$PHENOTYPE, levels = pheno_group_df$PHENOTYPE)
  df$VAR_TYPE <- factor(df$VAR_TYPE, levels = c("SNV", "INDEL", "DEL", "INS"))
  
  num_groups <- length(unique(pheno_group_df$GROUP))
  my_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "Pastel1"))(num_groups)
  
  group_bands <- pheno_group_df %>%
    group_by(GROUP) %>%
    summarize(ymin = min(PHENO_INDEX) - 0.5,
              ymax = max(PHENO_INDEX) + 0.5) %>%
    mutate(xmin = -Inf, xmax = Inf)
  
  plt <- ggplot(df, aes(x = POS, y = PHENO_INDEX)) +
      # Plot rectangles for phenotype groups
      geom_rect(data = group_bands, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = GROUP),
                  alpha = 0.2, color = NA, inherit.aes = FALSE) +
      geom_point(aes(shape = VAR_TYPE), size = 2, stroke = 0.6) +
      scale_x_continuous(
          labels = NULL,
          name = NULL,
          expand = expansion(mult = c(0.001, 0.001))
      ) +
      scale_y_continuous(
          breaks = pheno_group_df$PHENO_INDEX,
          labels = pheno_group_df$PHENOTYPE,
          name = "Phenotype"
      ) +
      scale_shape_manual(values = c(SNV = 16, INDEL = 16,    # filled circle (default)
                                      DEL = 15,      # filled square
                                      INS = 3,       # plus
                                      Other = 4))+
      scale_fill_manual(values = my_colors) +
  
      facet_grid(
          . ~ CHROM,
          scales = "free_x",
          space = "free_x",
          switch = "x"
      ) +
      coord_cartesian(clip = "off") + # Allow labels to extend outside the plot area
      theme_minimal() +
      theme(
          panel.grid = element_blank(),
          panel.spacing.x = unit(0.1, "lines"),
          strip.background = element_rect(fill = "gray95", color = "black"),
          strip.text = element_text(size = 9, face = "bold"),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "right",
          axis.title = element_blank(),
          plot.margin = margin(10, 10, 10, 10),
      ) +
          labs(
          shape = "Variant Type",
          fill = "Phenotype Group"
      )
  
  # Define type-specific output path
  type_plot_output_path <- sub("\\.pdf$", paste0("_", current_type, ".pdf"), plot_output_path)
  
  if (filter_to_sv) {
    type_rdata_output_path <- sub("\\.RData$", paste0("_", current_type, ".RData"), rdata_output_path)
    # Save the R plot object for further use
    save(plt, df, group_bands, file = type_rdata_output_path)
    cat(paste0("R plot object saved to: ", type_rdata_output_path, "\n"))
  } else {
    cat("No R plot object saved as filter_to_sv is FALSE.\n")
  }
  
  # Adjust height dynamically based on number of phenotypes
  num_phenos <- nrow(pheno_group_df)
  plot_height <- max(6, num_phenos * 0.25)
  
  pdf(type_plot_output_path, width = 15, height = plot_height)
  print(plt)
  dev.off()
  
  cat(paste0("Plot saved to: ", type_plot_output_path, "\n"))
}

