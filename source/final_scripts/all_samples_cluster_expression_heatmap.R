source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)

################################################################################

sce <- readRDS("data/cytof_objects/sclc_all_samples_with_clusters.rds")

all_markers <- readRDS("data/state_markers.rds")

markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam", "E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist", "SLUG", "PD-L1", "p-YAP", "CD44", "CD24")

ht <- create_expression_heatmap(sce, "new_clusters", markers_to_use,scale = T)

print(ht)



jpeg("figures/all_samples_cluster_expression_heatmap.jpg", width=200,height=100, units = "mm", res=1000)
print(ht)
dev.off()



