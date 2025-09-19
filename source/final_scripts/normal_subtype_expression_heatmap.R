################################################################################
# This script plots the protein expression from the "subtypes" found in normal 
# cells
################################################################################
source("source/sclc_cytof_functions.R")

set.seed(42)
################################################################################
# Read in data
################################################################################
ctcs <- readRDS("data/cytof_objects/normals_with_subtype.rds")


################################################################################
##Create heatmap
################################################################################
markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam", "E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist", "SLUG", "PD-L1", "p-YAP", "CD44", "CD24")

ht <- create_expression_heatmap(ctcs, "subtype", markers_to_use,"",scale = T)

################################################################################
# Save figure
################################################################################
tiff("figures/normal_subtype_expression_heatmap.tiff", width=200,height=100, units = "mm", res=1000)
print(ht)
dev.off()
