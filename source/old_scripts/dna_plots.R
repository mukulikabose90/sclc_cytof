library(EnvStats)
source("source/cytof_de_function.R")





sce <- readRDS("data/cytof_objects/sclc_cytof_sce_object_with_clusters_4.rds")


# Plot UMAP manually
xy <- reducedDim(sce, "UMAP")[, 1:2]
colnames(xy) <- c("x", "y")
df <- data.frame(colData(sce), xy, check.names = FALSE)






# Plot UMAP
ggplot(df)+
  geom_point(aes(x=x, y=y, color=new_clusters))+
  xlab("UMAP 1")+
  ylab("UMAP 2")+
  scale_color_manual(name = "Clusters",values=cluster_colors)+
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size=5)))+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=20), 
        axis.text = element_text(color = "black", size=15),
        axis.title = element_text(size=20),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18))



sce <- sce[,colData(sce)$experiment_id == 514521]


y <- assay(sce, "exprs")
df <- data.frame(t(y), colData(sce), check.names = FALSE)
gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(sce)))


colnames(gg_df)

temp <- gg_df %>% 
  dplyr::filter(antigen %in% c("Ir191","Ir193"))

ggplot(temp,aes(x = expression, color = condition)) + 
  facet_wrap(~antigen) + 
  geom_histogram(binwidth=.001)+
  ylab("normalized density") + 
  theme_classic() + 
  theme(panel.grid = element_blank(), 
                          strip.background = element_blank(), strip.text = element_text(face = "bold"), 
                          axis.text = element_text(color = "black"), axis.title = element_text(color = "black"))


ggplot(temp,aes(x = new_clusters, y= expression, color = new_clusters,)) + 
  geom_violin()
  facet_wrap(~antigen) + 
  ylab("normalized density") + 
  theme_classic() + 
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(), strip.text = element_text(face = "bold"), 
        axis.text = element_text(color = "black"), axis.title = element_text(color = "black"))





clusters <- sort(as.numeric(unique(colData(sce)$new_clusters)))
coefs <- c()
for(curr_cluster in clusters){
  
  curr_values <- gg_df %>% 
    dplyr::filter(new_clusters == curr_cluster & antigen == "Ir191") %>% 
    pull(expression)
  
  coefs <- append(coefs,cv(curr_values))
  
  
}

names(coefs) <- clusters

sort(coefs, decreasing = T)

boxplot(y["Ir191",colData(sce)$new_clusters == "normal"],y["Ir191",colData(sce)$condition == "cancer"])
wilcox.test(y["Ir191",colData(sce)$condition == "normal"],y["Ir191",colData(sce)$condition == "cancer"])


boxplot(y["Ir191",colData(sce)$condition == "normal"],y["Ir191",colData(sce)$condition == "cancer"])
wilcox.test(y["Ir191",colData(sce)$condition == "normal"],y["Ir191",colData(sce)$condition == "cancer"])

boxplot(y["Ir193",colData(sce)$condition == "normal"],y["Ir193",colData(sce)$condition == "cancer"])
wilcox.test(y["Ir193",colData(sce)$condition == "normal"],y["Ir193",colData(sce)$condition == "cancer"])

cv(y["Ir191",colData(sce)$condition == "normal"])
cv(y["Ir191",colData(sce)$condition == "cancer"])


cv(y["Ir193",colData(sce)$condition == "normal"])
cv(y["Ir193",colData(sce)$condition == "cancer"])







