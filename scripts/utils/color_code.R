# colors.R

library(colorspace)
library(ggsci)

# Base colors for each ancestry
base_colors <- c(
  "Holo" = "#D6C078",
  "Hakka" = "#99B581",
  "Southern Chinese" = "#4C87B8",
  "Northern Chinese" = "#9C5CBF",
  "Undefined Chinese" = "#C46A9C",
  "Other" = "#797D7F"
)

# Function to generate gradient colors
generate_gradient <- function(base_color) {
  colorRampPalette(c(base_color, lighten(base_color, 0.3)))(5)
}

# Expanded ancestry color definitions
ancestry_colors <- c(
  "Holo" = base_colors[["Holo"]],
  "Holo/Hakka" = generate_gradient(base_colors["Holo"])[2],
  "Holo/Southern Chinese" = generate_gradient(base_colors["Holo"])[3],
  "Holo/Northern Chinese" = generate_gradient(base_colors["Holo"])[4],
  "Holo/Undefined Chinese" = generate_gradient(base_colors["Holo"])[5],
  "Holo/Other" = generate_gradient(base_colors["Holo"])[5],
  
  "Hakka" = base_colors[["Hakka"]],
  "Hakka/Southern Chinese" = generate_gradient(base_colors["Hakka"])[2],
  "Hakka/Northern Chinese" = generate_gradient(base_colors["Hakka"])[3],
  "Hakka/Undefined Chinese" = generate_gradient(base_colors["Hakka"])[4],
  "Hakka/Other" = generate_gradient(base_colors["Hakka"])[5],
  
  "Southern Chinese" = base_colors[["Southern Chinese"]],
  "Southern Chinese/Northern Chinese" = generate_gradient(base_colors["Southern Chinese"])[2],
  "Southern Chinese/Undefined Chinese" = generate_gradient(base_colors["Southern Chinese"])[3],
  "Southern Chinese/Other" = generate_gradient(base_colors["Southern Chinese"])[4],
  
  "Northern Chinese" = base_colors[["Northern Chinese"]],
  "Northern Chinese/Undefined Chinese" = generate_gradient(base_colors["Northern Chinese"])[2],
  "Northern Chinese/Other" = generate_gradient(base_colors["Northern Chinese"])[3],
  
  "Undefined Chinese" = base_colors[["Undefined Chinese"]],
  "Other" = base_colors[["Other"]]
)

# Define the custom order
ancestry_custom_order <- c(
  "Holo", "Holo/Hakka", "Holo/Southern Chinese", "Holo/Northern Chinese", "Holo/Undefined Chinese", "Holo/Other",
  "Hakka", "Hakka/Southern Chinese", "Hakka/Northern Chinese", "Hakka/Undefined Chinese", "Hakka/Other",
  "Southern Chinese", "Southern Chinese/Northern Chinese", "Southern Chinese/Undefined Chinese", "Southern Chinese/Other",
  "Northern Chinese", "Northern Chinese/Undefined Chinese", "Northern Chinese/Other",
  "Undefined Chinese", "Other"
)

# Platform colors
platform_colors <- c("#1A759F", "#52B788", "#BE95C4")

################
# 1000GP colors
################
# Create population groups list
pop_groups <- list(
  AFR = c('MSL', 'ASW', 'LWK', 'ACB', 'ESN', 'GWD', 'YRI'),
  AMR = c('PUR', 'PEL', 'MXL', 'CLM'),
  EAS = c('JPT', 'KHV', 'CHB', 'CDX', 'CHS'),
  EUR = c('TSI', 'GBR', 'IBS', 'CEU', 'FIN'),
  SAS = c('BEB', 'STU', 'GIH', 'PJL', 'ITU')
)
super_pop <- names(pop_groups)
super_pop_colors <- setNames(pal_d3("category20")(20)[c(1:5)], super_pop)

pop <- unname(unlist(pop_groups))
# pop_colors <- setNames(c(pal_d3("category20")(20), pal_d3("category10")(6)), pop)

# Create population colors with variations of the superpopulation colors
# AFR populations - variations of blue
afr_colors <- c("#1A5276", "#377EB8", "#5499C7", "#7FB3D5", "#9ECAE1", "#C6DBEF", "#DEEBF7")
# AMR populations - variations of orange
amr_colors <- c("#FF7F00", "#FDB462", "#FED976", "#FEB24C", "#FED976")
# EAS populations - variations of green
eas_colors <- c("#196F3D", "#4DAF4A", "#74C476", "#A1D99B", "#C7E9C0")
# EUR populations - variations of red
eur_colors <- c("#E41A1C", "#FB6A4A", "#FC9272", "#FCBBA1", "#FEE0D2", "#FEE5D9", "#FCAE91")
# SAS populations - variations of purple
sas_colors <- c("#984EA3", "#AD7BB9", "#C2A5CF", "#D4B9DA", "#E7D4E8")

# Combine all colors in the order of populations
pop_colors <- setNames(
  c(afr_colors[1:length(pop_groups$AFR)],
    amr_colors[1:length(pop_groups$AMR)],
    eas_colors[1:length(pop_groups$EAS)],
    eur_colors[1:length(pop_groups$EUR)],
    sas_colors[1:length(pop_groups$SAS)]),
  pop
)
