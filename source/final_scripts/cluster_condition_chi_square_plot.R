source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################
sce <- readRDS("data/cytof_objects/sclc_all_samples_with_clusters.rds")

contin_table <- table(sce$new_clusters, sce$condition) %>% 
  as.data.frame.matrix()

test_res <- chisq.test(contin_table)

residuals <- test_res$stdres

# Calculate two-sided p-values from standard normal distribution
pvals <- 2 * (1 - pnorm(abs(residuals)))

# Adjust using Benjamini-Hochberg (FDR)
pvals_adj <- matrix(p.adjust(as.vector(pvals), method = "BH"),
                    nrow = nrow(pvals),
                    dimnames = dimnames(pvals))

# Convert to long format
res_df <- as.data.frame(as.table(residuals)) %>%
  mutate(p_value = as.vector(pvals_adj),
         stars = case_when(
           p_value < 0.001 ~ "***",
           p_value < 0.01  ~ "**",
           p_value < 0.05  ~ "*",
           TRUE            ~ ""))

colnames(res_df) <- c("cluster","condition","residual","adj_pval","stars")
res_df$condition <- ifelse(res_df$condition == "cancer", "Cancer","Normal")


res_df$condition <- factor(res_df$condition,levels = c("Normal","Cancer"))
res_df$cluster <- factor(res_df$cluster,levels = sort(levels(res_df$cluster),decreasing = T))

residual_plot <- ggplot(res_df, aes(x = condition, y = cluster)) +
  geom_point(aes(size = abs(residual), fill = residual), shape = 21, color = "black", stroke = 0.6) +
  geom_text(aes(label = stars), vjust = -1.5, size = 4) +  # Stars above dots
  scale_fill_gradient2(low = "#457B9D", mid = "white", high = "#dd4b33", midpoint = 0) +
  scale_size(range = c(2, 8)) +
  theme_minimal() +
  theme(axis.text = element_text(size=12,angle = 0, hjust = .5),
        axis.title = element_text(size=14)) +
  labs(title = "",
       size = "|Residual|", fill = "Residual",
       x = "",
       y = "")

residual_plot

################################################################################
# SAVE FIGURES
################################################################################

jpeg("figures/all_samples_cluster_condition_chi_results.jpg", width=100,height=180, units = "mm", res=1000)
print(residual_plot)
dev.off()


res_df$cluster <- factor(res_df$cluster,levels = sort(levels(res_df$cluster),decreasing = F))
cancer_enriched_clusters <- res_df %>% 
  filter(condition == "Cancer" & residual > 0) %>% 
  pull(cluster) %>% 
  unique() %>% 
  as.numeric()


saveRDS(cancer_enriched_clusters, "data/cancer_enriched_clusters.rds")







