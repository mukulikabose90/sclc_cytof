source("source/sclc_cytof_functions.R")
library(PCAtools)
set.seed(42)
################################################################################
# Read in data
################################################################################
sce <- readRDS("data/cytof_objects/sclc_all_samples_object_no_qc.rds")
sce <- sce[,sce$experiment_id != "531050"]

pca_result <- pca(sce@assays@data$exprs)

xy <- pca_result$rotated[,c(1,2)]

colnames(xy) <- c("pc1", "pc2")
df <- data.frame(colData(sce), xy, check.names = FALSE)

# Generate and save cluster colors
cluster_colors <- c(
  "#E57373",  # muted red
  "#FFB74D",  # muted orange
  "#90A4AE",  # slate gray
  "#81C784",  # muted green
  "#BA68C8",  # muted purple
  "#FFF176",  # soft yellow
  "#F48FB1",  # muted pink
  "#A1887F",  # warm tan
  "#64B5F6",  # muted blue
  "#4DB6AC",  # teal
  "#FFD54F",  # sunflower yellow
  "#CE93D8"   # lavender
)


# Plot PCA
p1 <- ggplot(df)+
  geom_point(aes(x=pc1, y=pc2, color=experiment_id),size=.1)+
  xlab("PC1")+
  ylab("PC2")+
  labs(color = "Experiment ID")+
  scale_color_manual(name = "Experiment ID",values=cluster_colors)+
  theme_classic() +
  guides(color = guide_legend(override.aes = list(size=5)))+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=8), 
        axis.text = element_text(color = "black", size=8),
        axis.title = element_text(size=8),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8))

p1

# Plot PCA
p2 <- ggplot(df)+
  geom_point(aes(x=pc1, y=pc2, color=experiment_id),size=.1)+
  facet_wrap(~experiment_id)+
  xlab("PC1")+
  ylab("PC2")+
  labs(color = "Experiment ID")+
  scale_color_manual(name = "Experiment ID",values=cluster_colors)+
  theme_classic() +
  guides(color = guide_legend(override.aes = list(size=5)))+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=8), 
        axis.text = element_text(color = "black", size=8),
        axis.title = element_text(size=8),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8))

p2


jpeg("figures/batch_pca.jpg", width=120,height=100, units = "mm", res=1000)
print(p1)
dev.off()

jpeg("figures/batch_pca_facet.jpg", width=200,height=120, units = "mm", res=1000)
print(p2)
dev.off()




