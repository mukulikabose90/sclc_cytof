
cluster_colors <- c("#dd4b33", "#F1FAEE", "#A8DADC", "#457B9D")
ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

sample_level_subtype <- ctcs@colData %>% 
  as.data.frame()  %>% 
  count(collection_id,subtype) %>% 
  group_by(collection_id) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n/total)*100) %>% 
  filter(total > 10) %>% 
  group_by(collection_id) %>% 
  filter(n == max(n)) %>% 
  select(collection_id,subtype) 


plot_df <- ctcs@colData %>% 
  as.data.frame() %>% 
  select(collection_id,treatment_status) %>% 
  merge(.,sample_level_subtype,by="collection_id") %>% 
  distinct() %>% 
  count(subtype,treatment_status) %>% 
  group_by(treatment_status) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n/total)*100)


plot_df$subtype <- factor(plot_df$subtype, levels = c("A","N","P","I"))

plot_df$treatment_status <- ifelse(plot_df$treatment_status == "naive", "Naive","Treated")
plot_df$treatment_status <- factor(plot_df$treatment_status, levels = c("Naive","Treated"))

ggplot(plot_df)+
  geom_col(aes(x=treatment_status,y=freq,fill=subtype))+
  geom_text(aes(x=treatment_status,y=102,label=total))+
  scale_fill_manual(values = cluster_colors)+
  theme_classic()+
  labs(x="",
       y="Percentage",
       fill="Subtype")
  
