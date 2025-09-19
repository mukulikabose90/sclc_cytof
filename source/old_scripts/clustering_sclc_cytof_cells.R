source("source/cytof_de_function.R")

script_seed <- 42
set.seed(script_seed)
################################################################################

sce <- readRDS("data/cytof_objects/sclc_cytof_sce_object.rds")

dim(sce)
# sce <- sce[,colData(sce)$sample_type != "cell_line"]
# 
# sce@metadata$experiment_info <- sce@metadata$experiment_info %>%
#   dplyr::filter(sample_type != "cell_line")

colData(sce)$condition <- factor(colData(sce)$condition, levels=c("normal", "cancer"))
sce@metadata$experiment_info$condition <- factor(sce@metadata$experiment_info$condition, levels=c("normal", "cancer"))

# sce <- sce[,!colData(sce)$patient_id %in% c("NORMAL16","NORMAL17")]


table(sce$condition)

sce <- CATALYST::cluster(sce, features = "state",
               xdim = 10, ydim = 10, maxK = 20, seed = script_seed)

# sce <- CATALYST::cluster(sce, features = "state",
#                          xdim = 10, ydim = 10, maxK = 20)

sce <- runDR(sce, "UMAP", cells = 5e3, features = "state")

sce@metadata$delta_area

plotDR(sce, "UMAP", color_by = "meta8", scale = T)+
  geom_point(size=1)+
  xlab("FFF")+
  scale_color_discrete(name = "New Legend Title")

plotDR(sce, "UMAP", color_by = "meta8", facet_by="condition",scale = T)

################################################################################

# Add new cluster assignments to colData
colData(sce)$new_clusters <- cluster_ids(sce, "meta8")

# Save data with cluster assignments
saveRDS(sce, "data/cytof_objects/sclc_cytof_sce_object_with_clusters.rds")

# Plot UMAP manually
xy <- reducedDim(sce, "UMAP")[, 1:2]
colnames(xy) <- c("x", "y")
df <- data.frame(colData(sce), xy, check.names = FALSE)



# Generate and save cluster colors
cluster_colors <- colorRampPalette(brewer.pal(10, "Paired"))(8)

display.brewer.pal(8,"Paired")

saveRDS(cluster_colors, "data/cluster_colors.rds")

cluster_colors <- readRDS("data/cluster_colors.rds")



# Plot UMAP
p1 <- ggplot(df)+
  geom_point(aes(x=x, y=y, color=new_clusters),size=.01)+
  xlab("UMAP 1")+
  ylab("UMAP 2")+
  scale_color_manual(name = "Clusters",values=cluster_colors)+
  theme_classic() +
  guides(color = guide_legend(override.aes = list(size=5)))+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=8), 
        axis.text = element_text(color = "black", size=8),
        axis.title = element_text(size=8),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8))



p1

facet_names <- c('normal'="Normal",'cancer'="Cancer")
p2 <- ggplot(df)+
  geom_point(aes(x=x, y=y, color=new_clusters), size=.01)+
  facet_wrap(~condition,labeller=as_labeller(facet_names))+
  xlab("UMAP 1")+
  ylab("UMAP 2")+
  scale_color_manual(name = "Clusters",values=cluster_colors)+
  theme_classic() +
  guides(color = guide_legend(override.aes = list(size=5)))+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=8), 
        axis.text = element_text(color = "black", size=8),
        axis.title = element_text(size=8),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8))


jpeg("figures/all_cells_cluster.jpg", width=120,height=100, units = "mm", res=1000)
print(p1)
dev.off()


jpeg("figures/all_cells_normal_vs_cancer_cluster.jpg", width=180,height=100, units = "mm", res=1000)
print(p2)
dev.off()





as.data.frame(colData(sce)) %>% 
  dplyr::count(condition)



