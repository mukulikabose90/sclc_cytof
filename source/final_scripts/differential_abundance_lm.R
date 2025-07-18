library(broom.mixed)

sce <- readRDS("data/cytof_objects/sclc_all_samples_with_clusters.rds")

data_df <- sce@colData %>% 
  as.data.frame() %>% 
  select(new_clusters,condition,sample_id) 

results_list <- list()
for(curr_cluster in sort(unique(curr_data$new_clusters))){
 
  data_df$cluster <- as.factor(as.integer(data_df$new_clusters == curr_cluster))
  
  formula_str <- glue("cluster ~ condition + (1 | sample_id)")
  
  
  model <- glmer(
    formula = as.formula(formula_str),
    family = binomial(link = "logit"),
    data = data_df)
  
  
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

summary(model)

plot_df <- all_results %>% 
  mutate(signif = ifelse(pval < 0.05, 16,1)) %>% 
  mutate(log_or = log(or)) %>% 
  mutate(log_upper_or = log(up_or)) %>% 
  mutate(log_lower_or = log(low_or)) %>% 
  mutate(star_height = ifelse(log_or > 0, log_or+.1,log_or-.1))



p1 <- ggplot(plot_df,aes(x=log_or,y=cluster))+
  geom_errorbarh(aes(xmin = log_lower_or, xmax = log_upper_or), height = 0.1, linewidth = .5)+
  geom_point(aes(x=log_or,y=cluster, shape = factor(signif),color=cluster),size=4,fill="white",show.legend = F)+
  scale_shape_manual(values = c("1" = 1, "16" = 16)) +
  geom_vline(xintercept = 0, linetype = 2)+
  xlim(-6,6)+
  scale_y_discrete(limits=rev)+
  labs(y="Cluster",
       x="log(OR)")+
  theme_classic()+
  theme(axis.text = element_text(size=12,angle = 0, hjust = 1),
        axis.title = element_text(size=14),
        axis.text.x = element_text(angle = 0, hjust = .5))


p1

jpeg("figures/all_samples_cluster_condition_or_results.jpg", width=100,height=180, units = "mm", res=1000)
print(p1)
dev.off()


cancer_enriched_clusters <- plot_df %>% 
  filter(log_or > 0 & pval < 0.05 ) %>% 
  pull(cluster) %>% 
  unique() %>% 
  as.numeric()


saveRDS(cancer_enriched_clusters, "data/cancer_enriched_clusters.rds")



