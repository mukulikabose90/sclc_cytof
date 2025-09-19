source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################
# Read in CyTOF data with cluster assignments
################################################################################
cancer_enriched <- readRDS("data/cytof_objects/cancer_enriched_with_clusters.rds")

#Subset to clusters with known CTC marker profiles
ctcs <- cancer_enriched[,cancer_enriched$new_clusters %in% c(2,3,4,5)]


ctcs %>% 
  colData() %>% 
  as.data.frame() %>% 
  count(new_clusters)

plot_df <- ctcs %>% 
  colData() %>% 
  as.data.frame() %>% 
  count(new_clusters, collection_id)


plot_df$new_clusters <- paste0("CTC Cluster ",plot_df$new_clusters)

ggplot(plot_df)+
  geom_col(aes(x=collection_id, y=n, fill = collection_id))+
  geom_text(aes(x=collection_id, y=n+5,label=n))+
  facet_wrap(~new_clusters, scales="free")+
  labs(x="Sample",
       y="Cell Count",
       fill="Sample")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))


