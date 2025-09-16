source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)

################################################################################

sce <- readRDS("data/cytof_objects/sclc_all_samples_with_clusters.rds")

all_markers <- readRDS("data/state_markers.rds")

markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam", "E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist", "SLUG", "PD-L1", "p-YAP", "CD44", "CD24")

scaled_heatmap <- create_expression_heatmap(sce, "new_clusters", markers_to_use,scale = T)

col_fun = colorRamp2(c(-3, -1, 0, 1, 3), 
                     c("#313695",  # deep blue
                       "#74add1",  # light blue
                       "#f7f7f7",  # white (center, 0)
                       "#f46d43",  # light red
                       "#a50026"))

ht <- Heatmap(t(scaled_heatmap),column_names_rot = 0,col = col_fun,
              cluster_columns = F, cluster_rows=F, column_title = "",
              row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 10),
              heatmap_legend_param = list(
                title = "   Scaled\nExpression",      
                title_gp = gpar(fontsize = 8), 
                labels_gp = gpar(fontsize = 6),
                legend_height = unit(3, "cm"),
                grid_width = unit(.5,"cm")))





print(ht)



tiff("figures/all_samples_cluster_expression_heatmap.tiff", width=200,height=60, units = "mm", res=600)
print(ht)
dev.off()



