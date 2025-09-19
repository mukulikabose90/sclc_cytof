source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
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

observed <- c()
expected <- c()

for(curr_cluster in clusters){
  
  a <- sum(colData(sce)$new_clusters == curr_cluster & colData(sce)$condition == "cancer")
  b <- sum(colData(sce)$new_clusters == curr_cluster & colData(sce)$condition != "cancer")
  c <- sum(colData(sce)$new_clusters != curr_cluster & colData(sce)$condition == "cancer")
  d <- sum(colData(sce)$new_clusters != curr_cluster & colData(sce)$condition != "cancer")
  
  contin_table <- matrix(c(a+.5,c+.5,b+.5,d+.5),ncol=2)
  
  chi_res <- chisq.test(contin_table)
  
  observed <- append(observed,a)
  expected <- append(expected,chi_res$expected[1,1])
  
  
  
  # ORs <- append(ORs,fisher_res$estimate)
  pvals <- append(pvals,chi_res$p.value)
}

fisher_res
(1697.20395*112.204)/(11.79605*16143.796)


pvals_adj <- p.adjust(pvals)

plot_df <- list()

plot_df[["value"]] <- c(observed,expected)
plot_df[["group"]] <- c(rep("observed", length(observed)),rep("expected", length(expected)))
plot_df[["cluster"]] <- c(clusters,clusters)
plot_df[["pval"]] <- c(pvals_adj,pvals_adj)

plot_df <- as.data.frame(plot_df)


plot_df <- plot_df %>% 
  mutate(significant = ifelse(pval < 0.05, "*","")) %>% 
  group_by(cluster) %>% 
  mutate(height = max(value))

condition_colors <- c("observed" = "#E57373","expected"="#64B5F6")

p1 <- ggplot(plot_df,aes(x=cluster,y=value,fill=group))+
  geom_col(position = "dodge")+
  geom_text(aes(y = height+.01,label=significant),size=6)+
  xlab("Cluster")+
  ylab("Count")+
  labs(fill="")+
  scale_fill_manual(values=condition_colors)+
  theme_classic()+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=8), 
        axis.text = element_text(color = "black", size=8),
        axis.title = element_text(size=8),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8))



p1

jpeg("figures/all_samples_cluster_diff_abundance_barplots.jpg", width=140,height=100, units = "mm", res=1000)
print(p1)
dev.off()








