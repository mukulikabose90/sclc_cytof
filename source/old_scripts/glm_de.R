
state_markers <- as.data.frame(rowData(sce)) %>% 
  dplyr::filter(marker_class == "state") %>% 
  pull(marker_name)

counts <- sce@assays@data$counts
counts <- counts[state_markers,]

clusters <- levels(colData(sce)$new_clusters)

protein <- c()
cluster <- c()
logfc <- c()
p_val <- c()

for(i in 1:length(state_markers)){
  
  # Set up data for fitting model
  data_i <- cbind(counts[i,],colData(sce))
  colnames(data_i)[1] <- "protein"
  
  # Fit glm
  fit <- glm("protein ~ new_clusters + experiment_id + patient_id + condition", data = data_i)

  # Extract p values for clusters
  cluster_pvals <- summary(fit)$coefficients[,4][1:length(clusters)]
  p_val <- append(p_val,cluster_pvals)
  
  for(curr_cluster in clusters){
    protein <- append(protein, state_markers[i])
    cluster <- append(cluster, curr_cluster)
    
    a <- counts[i,colData(sce)$new_clusters == curr_cluster]
    b <- counts[i,colData(sce)$new_clusters != curr_cluster]
    
    logfc <- append(logfc, log2(mean(a)/mean(b)))
    
  }
    
}


df <- data.frame(cbind(protein,cluster,logfc,p_val))

df$p_adj <- p.adjust(df$p_val, method = "BH")

df$logfc <- as.numeric(df$logfc)


plot_df <- df %>% 
  mutate(significance = factor(ifelse(p_adj < .05, "signif", "ns"), levels=c("signif","ns"))) %>% 
  dplyr::filter(cluster %in% signif_cluster) %>% 
  mutate(protein = reorder_within(protein, logfc, cluster))


ggplot(plot_df)+
  geom_col(aes(x=as.numeric(logfc), y=protein, fill=significance))+
  facet_wrap(~cluster, scales = "free_y")+
  scale_y_reordered()+
  ylab("Protein")+
  xlab("log2FC")









