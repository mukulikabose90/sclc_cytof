################################################################################
# This script fits a logistic regression mixed-effect model for each cluster to 
# determine which condition the cluster is enriched in (i.e. Cancer vs Normal).
# Then a forest plot displaying the log odds ratio is saved.
################################################################################
source("source/sclc_cytof_functions.R")

################################################################################
# Read in data
################################################################################
sce <- readRDS("data/cytof_objects/sclc_all_samples_with_clusters.rds")

################################################################################
# Fit logistic regression mixed-effect model for each cluster
################################################################################

data_df <- sce@colData %>% 
  as.data.frame() %>% 
  select(new_clusters,condition,sample_id) 

results_list <- list()
for(curr_cluster in sort(unique(data_df$new_clusters))){
 
  # binarize cluster 
  data_df$cluster <- as.factor(as.integer(data_df$new_clusters == curr_cluster))
  
  formula_str <- glue("cluster ~ condition + (1 | sample_id)")
  
  # fit model
  model <- glmer(
    formula = as.formula(formula_str),
    family = binomial(link = "logit"),
    data = data_df)
  
  # extract odds ratio, pvalue, CI and store in res
  cancer_or <- exp(fixef(model)["conditioncancer"])
  
  tidy_out <- tidy(model,effects='fixed')
  curr_pval <- tidy_out$p.value[2]
  
  lower_or <- exp(confint(model, parm = "beta_", method = "Wald")[2,][1])
  upper_or <- exp(confint(model, parm = "beta_", method = "Wald")[2,][2])
  
  res <- data.frame("cluster"=curr_cluster,"or"=cancer_or,"pval"=curr_pval,up_or=upper_or,low_or=lower_or)
  
  results_list <- append(results_list, list(res))
  
}

# Combine into one data frame
all_results <- bind_rows(results_list)

################################################################################
# Create plot dataframe
################################################################################

plot_df <- all_results %>% 
  mutate(padj = p.adjust(pval), method = "BH") %>% 
  mutate(signif = ifelse(padj < 0.05, 16,1)) %>% 
  mutate(log_or = log(or)) %>% 
  mutate(log_upper_or = log(up_or)) %>% 
  mutate(log_lower_or = log(low_or)) %>% 
  mutate(star_height = ifelse(log_or > 0, log_or+.1,log_or-.1))


p1 <- ggplot(plot_df,aes(x=log_or,y=fct_rev(cluster),color=cluster))+
  geom_point(aes(shape = factor(signif)),size=9,fill="white",show.legend = F)+
  scale_shape_manual(values = c("1" = 1, "16" = 16)) +
  geom_errorbarh(aes(xmin = log_lower_or, xmax = log_upper_or), height= 0.1, linewidth = .5,show.legend = F)+
  geom_vline(xintercept = 0, linetype = 2)+
  xlim(-5,5)+
  annotate("text", y=8.4, x=-2.5, label = "Normal", angle=0,size=5) +
  annotate("text", y=8.4, x=2.5, label = "Cancer", angle=0,size=5) +
  labs(y="Cluster",
       x="log(OR)")+
  theme_classic()+
  theme(axis.text = element_text(size=18,angle = 0, hjust = 1),
        axis.title = element_text(size=20),
        axis.text.x = element_text(angle = 0, hjust = .5))
p1
################################################################################
# Save figures
################################################################################
tiff("figures/all_samples_cluster_condition_or_results.tiff", width=100,height=200, units = "mm", res=600)
print(p1)
dev.off()

################################################################################
# Extract p-values
################################################################################
cluster_adj_pvalues <- as.numeric(sprintf("%.4f",plot_df$padj))
names(cluster_adj_pvalues) <- paste0("Cluster ", 1:nrow(plot_df))

cluster_adj_pvalues

################################################################################
# Extract cancer-enriched clusters
################################################################################
cancer_enriched_clusters <- plot_df %>% 
  filter(log_or > 0) %>% 
  pull(cluster) %>% 
  unique() %>% 
  as.numeric()

saveRDS(cancer_enriched_clusters, "data/cancer_enriched_clusters.rds")



