source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################

ctcs <- readRDS("data/cytof_objects/normals_with_subtype.rds")

table(ctcs$subtype)

markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam", "E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist", "SLUG", "PD-L1", "p-YAP", "CD44", "CD24")

ht <- create_expression_heatmap(ctcs, "subtype", markers_to_use,"",scale = T)

print(ht)

jpeg("figures/normal_subtype_expression_heatmap.jpg", width=200,height=100, units = "mm", res=1000)
print(ht)
dev.off()
