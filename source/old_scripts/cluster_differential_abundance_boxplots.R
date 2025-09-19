source("source/cytof_de_function.R")

script_seed <- 42
set.seed(script_seed)
################################################################################
# Differential abundance
FDR_cutoff <- 0.05

sce <- readRDS("data/cytof_objects/sclc_cytof_sce_object_with_clusters_4.rds")

cluster_prop_df <- as.data.frame(colData(sce)) %>%
  group_by(patient_id,new_clusters,condition) %>%
  summarise(n = n()) %>%
  mutate(total = sum(n)) %>% 
  group_by(patient_id) %>% 
  mutate(freq = (n / sum(n))*100)
  
cluster_prop_df$condition <- ifelse(cluster_prop_df$condition == "cancer", "Cancer", "Normal")
condition_colors <- c("Cancer" = "firebrick2","Normal"="royalblue")

p1 <- ggboxplot(cluster_prop_df, x="condition",y="freq",fill="condition", facet.by = "new_clusters",nrow=2)+
  stat_compare_means(size=4)+
  scale_fill_manual(values=condition_colors)+
  xlab("")+
  ylab("Percentage")+
  labs(fill="Condition")+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=15), 
        axis.text = element_text(color = "black", size=15),
        axis.title = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=20),
        legend.position="right")


p1

png("figures/cluster_diff_abundance_boxplots.png", width=16,height=8, units = "in", res=1200)
print(p1)
dev.off()



