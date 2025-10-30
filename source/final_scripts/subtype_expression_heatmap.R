################################################################################
# This script creates a heatmap displaying the scaled expression for each 
# protein in each of the 4 subtypes detected
################################################################################
source("source/sclc_cytof_functions.R")

set.seed(42)
################################################################################
# Read in data
################################################################################
ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

################################################################################
# Create protein x subtype heatmap 
################################################################################
markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam", "E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist", "SLUG", "PD-L1", "p-YAP", "CD44", "CD24")

scaled_heatmap <- create_expression_heatmap(ctcs, "subtype", markers_to_use,"", scale = T)


col_fun = colorRamp2(c(-3, -1, 0, 1, 3), 
                     c("#313695",  # deep blue
                       "#74add1",  # light blue
                       "#f7f7f7",  # white (center, 0)
                       "#f46d43",  # light red
                       "#a50026"))

ht <- Heatmap(t(scaled_heatmap),column_names_rot = 0,col = col_fun,
              cluster_columns = F, cluster_rows=F, column_title = "",column_names_centered = T,
              row_names_gp = gpar(fontsize = 16),column_names_gp = gpar(fontsize = 20),
              heatmap_legend_param = list(
                title = "   Scaled\nExpression",      
                title_gp = gpar(fontsize = 15), 
                labels_gp = gpar(fontsize = 14),
                legend_height = unit(3, "cm"),
                grid_width = unit(.5,"cm")))

################################################################################
# Save figure
################################################################################
tiff("figures/ctcs_subtype_expression_heatmap.tiff", width=200,height=100, units = "mm", res=600)
print(ht)
dev.off()

