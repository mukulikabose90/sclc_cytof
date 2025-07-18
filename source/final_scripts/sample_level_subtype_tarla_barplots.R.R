
cluster_colors <- c( "#F1FAEE", "#A8DADC", "#457B9D")
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
  select(collection_id,tarla) %>% 
  merge(.,sample_level_subtype,by="collection_id") %>% 
  distinct() %>% 
  count(subtype,tarla) %>% 
  group_by(tarla) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n/total)*100) %>% 
  filter(!is.na(tarla))


plot_df$subtype <- factor(plot_df$subtype, levels = c("A","N","P","I"))

plot_df$tarla <- ifelse(plot_df$tarla == "pre", "Pre-Tarlatamab","Post-Tarlatamab")
plot_df$tarla <- factor(plot_df$tarla, levels = c("Pre-Tarlatamab","Post-Tarlatamab"))

ggplot(plot_df)+
  geom_col(aes(x=tarla,y=freq,fill=subtype))+
  geom_text(aes(x=tarla,y=102,label=total))+
  scale_fill_manual(values = cluster_colors)+
  theme_classic()+
  labs(x="",
       y="Percentage",
       fill="Subtype")

