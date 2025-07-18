source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################
df <- readRDS("data/all_samples_umap_data.rds")

cancer_enriched_clusters <- readRDS("data/cancer_enriched_clusters.rds")

# plot with non-cancer_enricheds being translucent
df$cancer_enriched <- ifelse(df$new_clusters %in% cancer_enriched_clusters, "cancer_enriched", "non-cancer_enriched")

cluster_colors <- readRDS("data/cluster_colors.rds")

p1 <- ggplot(df)+
  geom_point(aes(x=x, y=y, color=new_clusters,alpha=cancer_enriched),size=.01)+
  xlab("UMAP 1")+
  ylab("UMAP 2")+
  labs(color = "Clusters")+
  # scale_color_manual(name = "Clusters",values=cluster_colors)+
  guides(color = guide_legend(override.aes = list(size=5)))+
  scale_alpha_manual(values = c("cancer_enriched" = 1, "non-cancer_enriched" = 0.05))+
  theme_classic() +
  guides(alpha = "none")+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=14), 
        strip.background = element_blank(),
        axis.text = element_text(color = "black", size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12))  


p1



facet_names <- c('normal'="Normal",'cancer'="Cancer")
p2 <- ggplot(df)+
  geom_point(aes(x=x, y=y, color=new_clusters,alpha=cancer_enriched),size=.01)+
  facet_wrap(~condition,labeller=as_labeller(facet_names))+
  xlab("UMAP 1")+
  ylab("UMAP 2")+
  labs(color = "Clusters")+
  # scale_color_manual(name = "Clusters",values=cluster_colors)+
  guides(color = guide_legend(override.aes = list(size=5)))+
  scale_alpha_manual(values = c("cancer_enriched" = 1, "non-cancer_enriched" = 0.05))+
  theme_classic() +
  guides(alpha = "none")+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=14), 
        strip.background = element_blank(),
        axis.text = element_text(color = "black", size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12))   


jpeg("figures/all_samples_cluster_opacity.jpg", width=120,height=100, units = "mm", res=1000)
print(p1)
dev.off()


jpeg("figures/all_samples_normal_vs_cancer_cluster_opacity.jpg", width=180,height=100, units = "mm", res=1000)
print(p2)
dev.off()