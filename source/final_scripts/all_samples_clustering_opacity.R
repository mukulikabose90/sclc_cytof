source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################

sce <- readRDS("data/cytof_objects/sclc_all_samples_object.rds")

ctc_clusters <- readRDS("data/ctc_clusters.rds")

dim(sce)

colData(sce)$condition <- factor(colData(sce)$condition, levels=c("normal", "cancer"))
sce@metadata$experiment_info$condition <- factor(sce@metadata$experiment_info$condition, levels=c("normal", "cancer"))

markers_to_use <- as.data.frame(rowData(sce)) %>% 
  dplyr::filter(marker_class == "state") %>% 
  pull(marker_name)


markers_to_use <- markers_to_use[-which(markers_to_use %in% c("NeuroD1","ASCL1","POU2F3","p-YAP","SLUG","Twist"))]

table(sce$condition)
dim(sce)

# sce <- CATALYST::cluster(sce, features = markers_to_use,
#                          xdim = 10, ydim = 10, maxK = 20, seed = script_seed)
# 
# sce <- runDR(sce, "UMAP", cells = 5e3, features = markers_to_use)

sce <- CATALYST::cluster(sce, features = "state",
                         xdim = 10, ydim = 10, maxK = 20, seed = script_seed)

sce <- runDR(sce, "UMAP", cells = 5e3, features = "state")

sce@metadata$delta_area

# plotDR(sce, color_by = "meta9",facet_by = "experiment_id")

################################################################################

# Add new cluster assignments to colData
colData(sce)$new_clusters <- cluster_ids(sce, "meta9")

# Plot UMAP manually
xy <- reducedDim(sce, "UMAP")[, 1:2]
colnames(xy) <- c("x", "y")
df <- data.frame(colData(sce), xy, check.names = FALSE)

# plot with non-ctcs being translucent
df$ctc <- ifelse(df$new_clusters %in% ctc_clusters, "ctc", "non-ctc")

cluster_colors <- readRDS("data/cluster_colors.rds")

p1 <- ggplot(df)+
  geom_point(aes(x=x, y=y, color=new_clusters,alpha=ctc),size=.01)+
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


p1



facet_names <- c('normal'="Normal",'cancer'="Cancer")
p2 <- ggplot(df)+
  geom_point(aes(x=x, y=y, color=new_clusters,alpha=ctc),size=.01)+
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


jpeg("figures/all_samples_cluster_opacity.jpg", width=120,height=100, units = "mm", res=1000)
print(p1)
dev.off()


jpeg("figures/all_samples_normal_vs_cancer_cluster_opacity.jpg", width=180,height=100, units = "mm", res=1000)
print(p2)
dev.off()