source("source/cytof_de_function.R")

script_seed <- 42
set.seed(script_seed)
################################################################################
# Differential abundance
FDR_cutoff <- 0.05

sce <- readRDS("data/cytof_objects/sclc_cytof_sce_object_with_clusters.rds")

ei <- metadata(sce)$experiment_info

da_formula2 <- createFormula(ei, 
                             cols_fixed = "condition", 
                             cols_random = c("sample_id", "patient_id", "experiment_id"))


contrast <- createContrast(c(0, 1))

da_res2 <- diffcyt(sce, 
                   formula = da_formula2, contrast = contrast,
                   analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",
                   clustering_to_use = "meta10", verbose = FALSE)





# Select significant clusters
signif_clusters <- data.frame(rowData(da_res2$res)) %>% 
  dplyr::filter(p_adj < FDR_cutoff) %>% 
  pull(cluster_id) %>% 
  as.numeric()


# Calculate proportion of each condition in each cluster
cluster_prop_df <- as.data.frame(colData(sce)) %>%
  group_by(new_clusters,condition) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))%>% 
  group_by(new_clusters) %>% 
  mutate(total = sum(n))

# Add star for significance
cluster_prop_df <- cluster_prop_df %>% 
  mutate(significant = ifelse(new_clusters %in% signif_clusters, "*",""))

# Plot abundance bar plots
ggplot(cluster_prop_df,aes(x=new_clusters,y=freq, fill=condition))+
  geom_col()+
  geom_text(aes(y = 1.03,label=total))+
  geom_text(aes(y = 1.05,label=significant),size=10)+
  xlab("Cluster")+
  ylab("Percentage")+
  theme()


ggplot(cluster_prop_df) +
  geom_boxplot(aes(x=new_clusters, y=freq))+
  facet_wrap(~condition)

