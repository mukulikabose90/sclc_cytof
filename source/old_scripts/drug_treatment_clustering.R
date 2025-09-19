source("source/cytof_de_function.R")

script_seed <- 42
set.seed(script_seed)
################################################################################

sce <- readRDS("data/cytof_objects/drug_treatment_experiment_sce_object.rds")

dim(sce)

sce <- CATALYST::cluster(sce, features = "state",
                         xdim = 10, ydim = 10, maxK = 20, seed = script_seed)


sce <- runDR(sce, "UMAP", cells = 5e3, features = "state")

sce@metadata$delta_area


plotDR(sce, "UMAP", color_by = "meta10", scale = T)

plotDR(sce, "UMAP", color_by = "meta10", facet_by = "day", scale = T)


# Add new cluster assignments to colData
colData(sce)$new_clusters <- cluster_ids(sce, "meta10")

################################################################################


clusters <- levels(colData(sce)$new_clusters)

pvals <- c()
ORs <- c()



for(curr_cluster in clusters){
  
  a <- sum(colData(sce)$new_clusters == curr_cluster & colData(sce)$day == "cancer")
  b <- sum(colData(sce)$new_clusters == curr_cluster & colData(sce)$day != "cancer")
  c <- sum(colData(sce)$new_clusters != curr_cluster & colData(sce)$day == "cancer")
  d <- sum(colData(sce)$new_clusters != curr_cluster & colData(sce)$day != "cancer")
  
  contin_table <- matrix(c(a+.5,c+.5,b+.5,d+.5),ncol=2)
  
  fisher_res <- fisher.test(contin_table)
  
  ORs <- append(ORs,fisher_res$estimate)
  pvals <- append(pvals,fisher_res$p.value)
}

# Select significant clusters
signif_clusters <- which(p.adjust(pvals) < 0.05)

cluster_prop_df <- as.data.frame(colData(sce)) %>% 
  dplyr::count(day,new_clusters) %>% 
  group_by(new_clusters) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n / total)*100) %>% 
  dplyr::filter(day !=0)

# Add star for significance
cluster_prop_df <- cluster_prop_df %>% 
  mutate(significant = ifelse(new_clusters %in% signif_clusters, "*","")) %>% 
  group_by(new_clusters) %>% 
  mutate(height = max(freq))


cluster_prop_df$day <- factor(cluster_prop_df$day, levels=c(0,3,7,10,14))

day_colors <- c("Cancer" = "firebrick2","Normal"="royalblue")

ggplot(cluster_prop_df,aes(x=new_clusters,y=freq,fill=day))+
  geom_col()+
  geom_text(aes(y = height+.01,label=significant),size=10)+
  xlab("Cluster")+
  ylab("Percentage")+
  labs(fill="day")+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=20), 
        axis.text = element_text(color = "black", size=15),
        axis.title = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=20))

################################################################################
cluster_prop_df <- as.data.frame(colData(sce)) %>% 
  dplyr::count(day,new_clusters) %>% 
  group_by(day) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n / total)*100) %>% 
  dplyr::filter(day !=0)

# Add star for significance
cluster_prop_df <- cluster_prop_df %>% 
  mutate(significant = ifelse(new_clusters %in% signif_clusters, "*","")) %>% 
  group_by(new_clusters) %>% 
  mutate(height = max(freq))


cluster_prop_df$day <- factor(cluster_prop_df$day, levels=c(0,3,7,10,14))

day_colors <- c("Cancer" = "firebrick2","Normal"="royalblue")

ggplot(cluster_prop_df,aes(x=day,y=freq,fill=new_clusters))+
  geom_col()+
  geom_text(aes(y = height+.01,label=significant),size=10)+
  xlab("Cluster")+
  ylab("Percentage")+
  labs(fill="day")+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=20), 
        axis.text = element_text(color = "black", size=15),
        axis.title = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=20))


