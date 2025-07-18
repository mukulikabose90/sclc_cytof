source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)

################################################################################
cancer_enriched <- readRDS("data/cytof_objects/cancer_enriched_with_clusters.rds")

table(cancer_enriched$new_clusters)

markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam", "E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist", "SLUG", "PD-L1", "p-YAP", "CD44", "CD24")

ht <- create_expression_heatmap(cancer_enriched, "new_clusters", markers_to_use)



jpeg("figures/cancer_enriched_cluster_expression_heatmap.jpg", width=200,height=100, units = "mm", res=1000)
print(ht)
dev.off()

