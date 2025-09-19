source("source/cytof_de_function.R")

script_seed <- 42
set.seed(script_seed)
################################################################################
# Differential abundance
FDR_cutoff <- 0.05

sce <- readRDS("data/cytof_objects/sclc_cytof_sce_object_with_clusters.rds")

sce <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")



sce <- sce[,colData(sce)$collection_id != "SC355-2"]
sce <- sce[,colData(sce)$condition == "cancer"]

sce <- ctcs

clusters <- levels(colData(sce)$new_clusters)

pvals <- c()
ORs <- c()

for(curr_cluster in clusters){
  
  a <- sum(colData(sce)$new_clusters == curr_cluster & colData(sce)$treatment_status == "cancer")
  b <- sum(colData(sce)$new_clusters == curr_cluster & colData(sce)$treatment_status != "cancer")
  c <- sum(colData(sce)$new_clusters != curr_cluster & colData(sce)$treatment_status == "cancer")
  d <- sum(colData(sce)$new_clusters != curr_cluster & colData(sce)$treatment_status != "cancer")
  
  contin_table <- matrix(c(a+.5,c+.5,b+.5,d+.5),ncol=2)
  
  fisher_res <- fisher.test(contin_table)
  
  ORs <- append(ORs,fisher_res$estimate)
  pvals <- append(pvals,fisher_res$p.value)
}

# Select significant clusters
signif_clusters <- which(p.adjust(pvals) < 0.05)

cluster_prop_df <- as.data.frame(colData(sce)) %>% 
  dplyr::count(new_clusters,treatment_status) %>% 
  group_by(treatment_status) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = n / total)

# Add star for significance
cluster_prop_df <- cluster_prop_df %>% 
  mutate(significant = ifelse(new_clusters %in% signif_clusters, "*","")) %>% 
  group_by(new_clusters) %>% 
  mutate(height = max(freq))


cluster_prop_df %>% 
  group_by(new_clusters) %>% 
  mutate(height = max(freq))


ggplot(cluster_prop_df,aes(x=new_clusters,y=freq,fill=treatment_status))+
  geom_col(position = "dodge")+
  geom_text(aes(y = height+.01,label=significant),size=10)+
  xlab("Cluster")+
  ylab("Percentage")




cluster_prop_df <- as.data.frame(colData(sce)) %>% 
  dplyr::count(subtype,treatment_status) %>% 
  group_by(treatment_status) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = n / total)

# Add star for significance
cluster_prop_df <- cluster_prop_df %>% 
  mutate(significant = ifelse(subtype %in% signif_clusters, "*","")) %>% 
  group_by(subtype) %>% 
  mutate(height = max(freq))


cluster_prop_df %>% 
  group_by(subtype) %>% 
  mutate(height = max(freq))


ggplot(cluster_prop_df,aes(x=treatment_status,y=freq,fill=subtype))+
  geom_col(position = "dodge")+
  geom_text(aes(y = height+.01,label=significant),size=10)+
  xlab("Cluster")+
  ylab("Percentage")





