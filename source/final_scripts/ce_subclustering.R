################################################################################
# This script re-clusters the cells from the initial cancer-enriched clusters 
# and plots the new clusters as a UMAP
################################################################################
source("source/sclc_cytof_functions.R")

set.seed(42)
################################################################################
# Read in CyTOF data with cluster assignments
################################################################################
sce <- readRDS("data/cytof_objects/sclc_all_samples_with_clusters.rds")

cancer_enriched_clusters <- readRDS("data/cancer_enriched_clusters.rds")

################################################################################
# Subset to cancer cells in cancer_enriched cluster
################################################################################
cancer_enriched <- sce[,colData(sce)$new_clusters %in% cancer_enriched_clusters]
cancer_enriched <- cancer_enriched[,colData(cancer_enriched)$condition == "cancer"]

################################################################################
# Cluster cells using FlowSOM then run UMAP
################################################################################
markers_to_use <- readRDS("data/state_markers.rds")
# markers_to_use <- markers_to_use[-which(markers_to_use %in% c("NeuroD1","ASCL1","POU2F3","p-YAP","SLUG","Twist","p-Rb"))]

cancer_enriched <- CATALYST::cluster(cancer_enriched, features = markers_to_use,
                          xdim = 10, ydim = 10, maxK = 20, seed = 42)

cancer_enriched <- runDR(cancer_enriched, "UMAP", cells = 5e3, features = markers_to_use)

################################################################################
# Identify optimal number of clusters and assign cells
################################################################################
cancer_enriched@metadata$delta_area

CATALYST::plotDR(cancer_enriched, color_by = "meta5",facet_by = "condition")
CATALYST::plotDR(cancer_enriched, color_by = "meta6")

colData(cancer_enriched)$new_clusters <- cluster_ids(cancer_enriched, "meta8")

# Save data with cluster assignments
saveRDS(cancer_enriched, "data/cytof_objects/cancer_enriched_with_clusters.rds")

################################################################################
# Create UMAP figures
################################################################################
xy <- reducedDim(cancer_enriched, "UMAP")[, 1:2]
colnames(xy) <- c("x", "y")
df <- data.frame(colData(cancer_enriched), xy, check.names = FALSE)

# Generate and save cluster colors
cluster_colors <- c(
  "#E57373",  # muted red
  "#FFB74D",  # muted orange
  "#81C784",  # muted green
  "#BA68C8",  # muted purple
  "#FFF176",  # soft yellow
  "#F48FB1",  # muted pink
  "#A1887F",  # warm tan
  "#64B5F6",  # muted blue
  "#4DB6AC"   # soft teal
)

saveRDS(cluster_colors, "data/cluster_colors.rds")

cluster_colors <- readRDS("data/cluster_colors.rds")

# Plot UMAP
p1 <- ggplot(df)+
  geom_point(aes(x=x, y=y, color=new_clusters),size=.1)+
  xlab("UMAP 1")+
  ylab("UMAP 2")+
  labs(color = "Clusters")+
  # scale_color_manual(name = "Clusters",values=cluster_colors)+
  theme_classic() +
  guides(color = guide_legend(override.aes = list(size=5)))+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=8), 
        axis.text = element_text(color = "black", size=8),
        axis.title = element_text(size=8),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8))

################################################################################
# Save figures
################################################################################
tiff("figures/cancer_enriched_cluster.tiff", width=120,height=100, units = "mm", res=1000)
print(p1)
dev.off()


