source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################

treatment_status_sce <- readRDS("data/cytof_objects/treatment_status_sce.rds")

markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam", "E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist", "SLUG", "PD-L1", "p-YAP", "CD44", "CD24")

ht <- create_expression_heatmap(treatment_status_sce, "new_clusters", markers_to_use)

jpeg("figures/treatment_status_cluster_expression_heatmap.jpg", width=300,height=100, units = "mm", res=1000)
print(ht)
dev.off()
