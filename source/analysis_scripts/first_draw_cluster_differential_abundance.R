source("source/cytof_de_function.R")

script_seed <- 42
set.seed(script_seed)
################################################################################
# Differential abundance
FDR_cutoff <- 0.05

sce <- readRDS("data/cytof_objects/sclc_first_draw_with_clusters.rds")

clusters <- levels(colData(sce)$new_clusters)

pvals <- c()
ORs <- c()

for(curr_cluster in clusters){
  
  a <- sum(colData(sce)$new_clusters == curr_cluster & colData(sce)$condition == "cancer")
  b <- sum(colData(sce)$new_clusters == curr_cluster & colData(sce)$condition != "cancer")
  c <- sum(colData(sce)$new_clusters != curr_cluster & colData(sce)$condition == "cancer")
  d <- sum(colData(sce)$new_clusters != curr_cluster & colData(sce)$condition != "cancer")
  
  contin_table <- matrix(c(a+.5,c+.5,b+.5,d+.5),ncol=2)
  
  fisher_res <- fisher.test(contin_table)
  
  ORs <- append(ORs,fisher_res$estimate)
  pvals <- append(pvals,fisher_res$p.value)
}

# Select significant clusters
signif_clusters <- which(p.adjust(pvals) < 0.05)

cluster_prop_df <- as.data.frame(colData(sce)) %>% 
  dplyr::count(new_clusters,condition) %>% 
  group_by(condition) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n / total)*100)

# Add star for significance
cluster_prop_df <- cluster_prop_df %>% 
  mutate(significant = ifelse(new_clusters %in% signif_clusters, "*","")) %>% 
  group_by(new_clusters) %>% 
  mutate(height = max(freq))


cluster_prop_df$condition <- ifelse(cluster_prop_df$condition == "cancer", "Cancer", "Normal")

condition_colors <- c("Cancer" = "firebrick2","Normal"="royalblue")

p1 <- ggplot(cluster_prop_df,aes(x=new_clusters,y=freq,fill=condition))+
  geom_col(position = "dodge")+
  geom_text(aes(y = height+.01,label=significant),size=6)+
  xlab("Cluster")+
  ylab("Percentage")+
  labs(fill="Condition")+
  scale_fill_manual(values=condition_colors)+
  theme_classic()+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=8), 
        axis.text = element_text(color = "black", size=8),
        axis.title = element_text(size=8),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8))




p1

jpeg("figures/first_draw_cluster_diff_abundance_barplots.jpg", width=120,height=100, units = "mm", res=1000)
print(p1)
dev.off()

