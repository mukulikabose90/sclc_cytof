source("source/sclc_cytof_functions.R")

set.seed(42)
FDR_cutoff <- 0.05
################################################################################
# Read in data
################################################################################

sce <- readRDS("data/cytof_objects/sclc_all_samples_with_clusters.rds")

################################################################################
# Calculate differential abundance using Fisher's exact test
################################################################################
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

################################################################################
# Create DA barplot 
################################################################################

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

# reorder condition factor
cluster_prop_df$condition <- ifelse(cluster_prop_df$condition == "cancer", "Cancer", "Normal")
cluster_prop_df$condition <- factor(cluster_prop_df$condition, levels=c("Normal","Cancer"))
  
condition_colors <- c("Cancer" = "#E57373","Normal"="#64B5F6")

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

jpeg("figures/all_samples_cluster_diff_abundance_barplots.jpg", width=120,height=100, units = "mm", res=1000)
print(p1)
dev.off()

################################################################################
# Identify and save CTC clusters
################################################################################
ctc_clusters <- c()
for(i in unique(cluster_prop_df$new_clusters)){
  cancer_freq <- cluster_prop_df %>% 
    dplyr::filter(significant == "*" & new_clusters == i & condition == "Cancer") %>%
    pull(freq)
  
  normal_freq <- cluster_prop_df %>% 
    dplyr::filter(significant == "*" & new_clusters == i & condition == "Normal") %>%
    pull(freq)
  
  if(length(cancer_freq) == 0){
    cancer_freq <- 0
  }
  
  if(length(normal_freq) == 0){
    normal_freq <- 0
  }
  
  if(cancer_freq > normal_freq){
    ctc_clusters <- append(ctc_clusters, i)
  }
}

saveRDS(as.numeric(ctc_clusters), "data/ctc_clusters.rds")

