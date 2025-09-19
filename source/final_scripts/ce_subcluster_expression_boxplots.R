################################################################################
# This script generates boxplots for selected proteins displaying the 
# distribution of expression across all cells in the given comparison groups
################################################################################
source("source/sclc_cytof_functions.R")

set.seed(42)
################################################################################
# Read in data
################################################################################
cancer_enriched <- readRDS("data/cytof_objects/cancer_enriched_with_clusters.rds")

################################################################################
# Plot boxplots for SCLC TFs expression from clusters
################################################################################
markers_to_use <- c("NeuroD1","ASCL1","POU2F3")

cancer_enriched$ctc <- ifelse(cancer_enriched$new_clusters %in% c(1,3),"Non-CTCs", "CTCs")

p1 <- create_marker_boxplots(cancer_enriched,markers_to_use,"new_clusters",fill = "ctc", nrow_to_use=1)

p1 <- p1+
  labs(y="Expression",
       x="",
       fill="")+
  scale_fill_manual(
    values = c("CTCs" = "#f46d43","Non-CTCs" = "#74add1"))+
  rremove("legend")

################################################################################
# Plot boxplots for selected proteins expression from clusters
################################################################################
# markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam", "E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist", "SLUG", "PD-L1", "p-YAP", "CD44", "CD24")
markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam", "E-Cad", "MUC-1", "Vimentin", "PD-L1", "CD24")

p2 <- create_marker_boxplots(cancer_enriched,markers_to_use,"new_clusters",fill="ctc", nrow_to_use=2)

p2 <- p2+
  labs(y="Expression",
       x="",
       fill="")+
  scale_fill_manual(
    values = c("CTCs" = "#f46d43","Non-CTCs" = "#74add1"))+
  rremove("legend")

################################################################################
# Save figures
################################################################################
tiff(glue("figures/cancer_enriched_cluster_tf_expression_boxplot.tiff"), width=160,height=160, units = "mm", res=1000)
print(p1)
dev.off()

tiff(glue("figures/cancer_enriched_cluster_expression_boxplot.tiff"), width=300,height=200, units = "mm", res=1000)
print(p2)
dev.off()

