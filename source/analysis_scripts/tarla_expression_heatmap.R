source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################

tarla_sce <- readRDS("data/cytof_objects/tarla_sce.rds")

markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam", "E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist", "SLUG", "PD-L1", "p-YAP", "CD44", "CD24")

ht <- create_expression_heatmap(tarla_sce, "new_clusters", all_markers)

jpeg("figures/tarla_cluster_expression_heatmap.jpg", width=300,height=100, units = "mm", res=1000)
print(ht)
dev.off()
