source("source/sclc_cytof_functions.R")

set.seed(42)
################################################################################
# Read in data
################################################################################
sce <- readRDS("data/cytof_objects/sclc_all_samples_object.rds")

# Reorder factors
colData(sce)$condition <- factor(colData(sce)$condition, levels=c("normal", "cancer"))
sce@metadata$experiment_info$condition <- factor(sce@metadata$experiment_info$condition, levels=c("normal", "cancer"))

################################################################################
# Get and save state markers
################################################################################
marker_info <- read.csv("data/cytof_panel_info.csv")
marker_info <- data.frame(marker_info, stringsAsFactors = FALSE)

state_markers <- marker_info %>%
  dplyr::filter(marker_class == "state") %>%
  pull(antigen)

saveRDS(state_markers, "data/state_markers.rds")
################################################################################
# Cluster cells using FlowSOM then run UMAP
################################################################################
markers_to_use <- state_markers
# markers_to_use <- markers_to_use[-which(markers_to_use %in% c("NeuroD1","ASCL1","POU2F3","p-YAP","SLUG","Twist"))]

sce <- CATALYST::cluster(sce, features = markers_to_use,
                         xdim = 10, ydim = 10, maxK = 20, seed = script_seed)

sce <- runDR(sce, "UMAP", cells = 5e3, features = markers_to_use)

# sce <- CATALYST::cluster(sce, features = "state",
#                          xdim = 10, ydim = 10, maxK = 20, seed = script_seed)
# 
# sce <- runDR(sce, "UMAP", cells = 5e3, features = "state")

################################################################################
# Identify optimal number of clusters and assign cells
################################################################################
sce@metadata$delta_area

CATALYST::plotDR(sce, color_by = "meta9",facet_by = "experiment_id")

CATALYST::plotDR(sce, color_by = "experiment_id")

# Add new cluster assignments to colData
colData(sce)$new_clusters <- cluster_ids(sce, "meta9")

# Save data with cluster assignments
saveRDS(sce, "data/cytof_objects/sclc_all_samples_with_clusters.rds")

################################################################################
# Create UMAP figures
################################################################################
xy <- reducedDim(sce, "UMAP")[, 1:2]
colnames(xy) <- c("x", "y")
df <- data.frame(colData(sce), xy, check.names = FALSE)

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
  
p1

# Plot UMAP faceted by condition
facet_names <- c('normal'="Normal",'cancer'="Cancer")
p2 <- ggplot(df)+
  geom_point(aes(x=x, y=y, color=new_clusters),size=.01)+
  facet_wrap(~condition,labeller=as_labeller(facet_names))+
  xlab("UMAP 1")+
  ylab("UMAP 2")+
  labs(color = "Clusters")+
  # scale_color_manual(name = "Clusters",values=cluster_colors)+
  guides(color = guide_legend(override.aes = list(size=5)))+
  scale_alpha_manual(values = c("ctc" = 1, "non-ctc" = 0.05))+
  theme_classic() +
  guides(alpha = "none")+
  theme(panel.grid.minor = element_blank(), 
      strip.text = element_text(face = "bold", size=8), 
      axis.text = element_text(color = "black", size=8),
      axis.title = element_text(size=8),
      legend.text = element_text(size=6),
      legend.title = element_text(size=8))  

p2

################################################################################
# Save figures
################################################################################
jpeg("figures/all_samples_cluster.jpg", width=120,height=100, units = "mm", res=1000)
print(p1)
dev.off()


jpeg("figures/all_samples_normal_vs_cancer_cluster.jpg", width=180,height=100, units = "mm", res=1000)
print(p2)
dev.off()
