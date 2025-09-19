source("source/cytof_de_function.R")

script_seed <- 42
set.seed(script_seed)
################################################################################

sce <- readRDS("data/cytof_objects/sclc_cytof_sce_object_with_clusters.rds")

cluster_prop_df <- as.data.frame(colData(sce)) %>%
  group_by(patient_id,new_clusters) %>%
  summarise(n = n()) %>%
  mutate(freq = (n / sum(n))*100)%>% 
  group_by(patient_id) %>% 
  mutate(total = sum(n))


cluster_prop_df$condition <- ifelse(grepl("SC", cluster_prop_df$patient_id), "Cancer","Normal")

# read in cluster colors
cluster_colors <- readRDS("data/cluster_colors.rds")

# Set opacity of cluster to .9, except "CTC" clusters
cluster_prop_df$alpha <- ifelse(cluster_prop_df$new_clusters %in% c(1,2,5,6),1,0.99)

cluster_prop_df <- cluster_prop_df %>% 
  dplyr::filter(total > 29)


p1 <- ggplot(cluster_prop_df,aes(x=patient_id,y=freq, fill=new_clusters))+
  geom_col()+
  geom_text(aes(y = 103,label=total))+
  facet_grid(~condition, scales = "free_x",space='free')+
  scale_fill_manual(name = "Clusters",values=cluster_colors)+
  scale_alpha(guide = 'none')+
  xlab("Patient")+
  ylab("Percentage")+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=20), 
        axis.text = element_text(color = "black", size=15),
        axis.text.x = element_text(angle=70,hjust=1,vjust=1),
        axis.title = element_text(size=20),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18))

p2 <- ggplot(cluster_prop_df,aes(x=patient_id,y=freq, fill=new_clusters, alpha=alpha))+
  geom_col()+
  geom_text(aes(y = 103,label=total))+
  facet_grid(~condition, scales = "free_x",space='free')+
  scale_fill_manual(name = "Clusters",values=cluster_colors)+
  scale_alpha(guide = 'none')+
  xlab("Patient")+
  ylab("Percentage")+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=20), 
        axis.text = element_text(color = "black", size=15),
        axis.text.x = element_text(angle=70,hjust=1,vjust=1),
        axis.title = element_text(size=20),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18))
  

p1

p2

# png("figures/patient_cluster_proportions_all.png", width=16,height=8, units = "in", res=1200)
# print(p1)
# dev.off()

png("figures/patient_cluster_proportions_signif_without_N16_N17.png", width=16,height=8, units = "in", res=1200)
print(p2)
dev.off()


################################################################################
# Compare distribution of proportions across patients

p3 <- ggboxplot(cluster_prop_df %>% 
            dplyr::filter(new_clusters == 1), x="condition",
          y= "freq",fill="condition",lwd=1)+
  scale_fill_manual(values=c("white", "white"))+
  xlab("Condition")+
  ylab("Percentage")+
  stat_compare_means(size=8, label.y=.6, label.x=.75)+
  theme(legend.position="none",
        axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(size=1))


p3
# png("figures/cluster_1_patient_proportions_no_color.png", width = 10, height = 10, units = "in", res = 1000)
tiff("figures/cluster_1_patient_proportions_no_color.tiff", width = 10, height = 10, units = "in", res = 1000)
print(p3)
dev.off()

# ################################################################################
# # Check proportions in each experiment
# cluster_prop_df <- as.data.frame(colData(sce))
# 
# cluster_prop_df <- cluster_prop_df %>%
#   group_by(experiment_id,new_clusters) %>%
#   summarise(n = n()) %>%
#   mutate(freq = n / sum(n))%>% 
#   group_by(experiment_id) %>% 
#   mutate(total = sum(n))
# 
# 
# ggplot(cluster_prop_df,aes(x=experiment_id,y=freq, fill=new_clusters))+
#   geom_col()+
#   geom_text(aes(y = 1.03,label=total))+
#   xlab("Experiment")+
#   ylab("Percentage")+
#   theme(axis.text.x = element_text(angle=70,hjust=1,vjust=1))
# 
# 
# ggplot(cluster_prop_df,aes(x=experiment_id,y=freq, fill=new_clusters))+
#   geom_col()+
#   geom_text(aes(y = 1.03,label=total))+
#   xlab("Experiment")+
#   ylab("Percentage")+
#   theme(axis.text.x = element_text(angle=70,hjust=1,vjust=1))+
#   scale_fill_manual(values=ifelse(c(1:10) == 1, "blue","gray"))
# 
# ################################################################################
# # Check proportions in each sample_id
# cluster_prop_df <- as.data.frame(colData(sce))
# 
# cluster_prop_df <- cluster_prop_df %>%
#   group_by(sample_id,new_clusters) %>%
#   summarise(n = n()) %>%
#   mutate(freq = n / sum(n))%>% 
#   group_by(sample_id) %>% 
#   mutate(total = sum(n))
# 
# cluster_prop_df$condition <- ifelse(grepl("SC", cluster_prop_df$sample_id), "cancer","normal")
# 
# ggplot(cluster_prop_df,aes(x=sample_id,y=freq, fill=new_clusters))+
#   geom_col()+
#   geom_text(aes(y = 1.03,label=total))+
#   xlab("Experiment")+
#   ylab("Percentage")+
#   theme(axis.text.x = element_text(angle=70,hjust=1,vjust=1))
# 
# 
# ggplot(cluster_prop_df,aes(x=sample_id,y=freq, fill=new_clusters))+
#   geom_col()+
#   geom_text(aes(y = 1.03,label=total))+
#   facet_grid(~condition, scales = "free_x")+
#   xlab("Sample")+
#   ylab("Percentage")+
#   theme(axis.text.x = element_text(angle=70,hjust=1,vjust=1))+
#   scale_fill_manual(values=ifelse(c(1:10) == 1, "blue","gray"))






