
# Save data with subtype assignments
ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

################################################################################
# Plot proportion of subtypes in each collection 
cluster_prop_df <- as.data.frame(colData(ctcs)) %>%
  group_by(collection_id,subtype) %>%
  summarise(n = n()) %>%
  mutate(freq = (n / sum(n))*100)%>% 
  group_by(collection_id) %>% 
  mutate(total = sum(n))

cluster_prop_df <- cluster_prop_df %>% 
  dplyr::filter(total > 29)


subtype_colors <- colorRampPalette(brewer.pal(4,"Set1"))(4)

p1 <- ggplot(cluster_prop_df,aes(x=collection_id, y=freq, fill=subtype))+
  geom_col()+
  geom_text(aes(y = 103,label=total))+
  xlab("Sample")+
  ylab("Percentage")+
  labs(fill="Subtype")+
  scale_fill_manual(values= subtype_colors)+
  theme_classic()+
  theme(axis.text.x = element_text(angle=70,hjust=1,vjust=1))+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=20), 
        axis.text = element_text(color = "black", size=15),
        axis.title = element_text(size=20),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18))

png("figures/ctcs_subtype_proportions.png", width = 12, height = 6, units = "in", res = 1200)
print(p1)
dev.off()

################################################################################
# Plot proportion of subtypes in each time point 
cluster_prop_df <- as.data.frame(colData(ctcs)) %>%
  group_by(collection_id,subtype,patient_id,treatment_status,sample_num) %>%
  summarise(n = n()) %>%
  group_by(collection_id) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = n/total)

cluster_prop_df <- cluster_prop_df %>% 
  dplyr::filter(total > 5)

ggplot(cluster_prop_df,aes(x=treatment_status, y=freq, fill=subtype))+
  facet_wrap(~patient_id, scales='free_x')+
  geom_col(position="dodge")+
  geom_text(aes(y = 1.03,label=total))+
  xlab("Sample")+
  ylab("Percentage")+
  theme(axis.text.x = element_text(angle=70,hjust=1,vjust=1))

################################################################################
# Treatment status proportions in each subtype
clusters <- unique(colData(ctcs)$subtype)

pvals <- c()
ORs <- c()

curr_cluster <- "N"
198.5/520.5

1313.5/4266.5
for(curr_cluster in clusters){
  
  a <- sum(colData(ctcs)$subtype == curr_cluster & colData(ctcs)$treatment_status == "naive")
  b <- sum(colData(ctcs)$subtype == curr_cluster & colData(ctcs)$treatment_status != "naive")
  c <- sum(colData(ctcs)$subtype != curr_cluster & colData(ctcs)$treatment_status == "naive")
  d <- sum(colData(ctcs)$subtype != curr_cluster & colData(ctcs)$treatment_status != "naive")
  
  contin_table <- matrix(c(a+.5,c+.5,b+.5,d+.5),ncol=2)
  
  fisher_res <- fisher.test(contin_table)
  
  ORs <- append(ORs,fisher_res$estimate)
  pvals <- append(pvals,fisher_res$p.value)
}

# Select significant clusters
signif_clusters <- clusters[which(p.adjust(pvals) < 0.05)]

cluster_prop_df <- as.data.frame(colData(ctcs)) %>% 
  dplyr::filter(treatment_status %in% c("naive","treated")) %>% 
  dplyr::count(subtype,treatment_status) %>% 
  group_by(treatment_status) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n / total)*100)

# Add star for significance
cluster_prop_df <- cluster_prop_df %>% 
  mutate(significant = ifelse(subtype %in% signif_clusters, "*","")) %>% 
  group_by(subtype) %>% 
  mutate(height = max(freq))

# Capitalize status
cluster_prop_df$treatment_status <- ifelse(cluster_prop_df$treatment_status == "naive", "Naive","Treated")

p2 <- ggplot(cluster_prop_df,aes(x=subtype,y=freq,fill=treatment_status))+
  geom_col(position = "dodge")+
  geom_text(aes(y = height+.01,label=significant),size=10)+
  scale_fill_manual(values=c("deepskyblue3", "darkorange"))+
  labs(fill="Treatment Status")+
  xlab("Subtype")+
  ylab("Percentage")+
  theme_classic()+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=20), 
        axis.text = element_text(color = "black", size=15),
        axis.title = element_text(size=20),
        axis.ticks = element_line(linewidth =1),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18))

png("figures/treatment_status_subtype_barplots.png", width = 7, height = 6, units = "in", res = 1200)
print(p2)
dev.off()





