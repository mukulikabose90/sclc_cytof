source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

cluster_colors <- c("#dd4b33", "#F1FAEE", "#A8DADC", "#457B9D")

################################################################################
# Plot total subtype proportions between naive and treated
################################################################################
curr_data <- ctcs

plot_df <- as.data.frame(curr_data@colData) %>% 
  filter(tarla != "post" | is.na(tarla)) %>% 
  select(collection_id,treatment_status,subtype) %>% 
  dplyr::count(collection_id,treatment_status,subtype) %>% 
  group_by(collection_id,treatment_status) %>% 
  mutate(total = sum(n)) %>% 
  filter(total > 10) %>% 
  mutate(freq = (n/total)*100) 

plot_df$subtype <- factor(plot_df$subtype, levels=c("A","N","P",'I'))

plot_df$treatment_status <- ifelse(plot_df$treatment_status == "naive","Naive","Treated")

ggplot(plot_df,aes(x=subtype,y=freq))+
  geom_violin(aes(fill=subtype),draw_quantiles = .5,)+
  geom_jitter()+
  facet_wrap(~treatment_status)+
  scale_fill_manual(name = "Subtype",values=cluster_colors)+
  theme_classic()


ggplot(plot_df,aes(x=treatment_status,y=freq))+
  geom_violin(aes(fill=subtype),draw_quantiles = .5,)+
  stat_compare_means()+
  geom_jitter()+
  facet_wrap(~subtype)+
  scale_fill_manual(name = "Subtype",values=cluster_colors)+
  theme_classic()



################################################################################

curr_data <- ctcs

plot_df <- as.data.frame(curr_data@colData) %>% 
  filter(!is.na(tarla)) %>% 
  select(collection_id,tarla,subtype) %>% 
  dplyr::count(collection_id,tarla,subtype) %>% 
  group_by(collection_id,tarla) %>% 
  mutate(total = sum(n)) %>% 
  filter(total > 10) %>% 
  mutate(freq = (n/total)*100) 

plot_df$subtype <- factor(plot_df$subtype, levels=c("A","N","P",'I'))

plot_df$tarla <- ifelse(plot_df$tarla == "pre","Pre-Tarlatamab","Post-Tarlatamab")
plot_df$tarla <- factor(plot_df$tarla, levels=c("Pre-Tarlatamab","Post-Tarlatamab"))

ggplot(plot_df,aes(x=subtype,y=freq))+
  geom_violin(aes(fill=subtype),draw_quantiles = .5,)+
  geom_jitter()+
  facet_wrap(~tarla)+
  scale_fill_manual(name = "Subtype",values=cluster_colors)+
  theme_classic()


ggplot(plot_df,aes(x=tarla,y=freq))+
  geom_violin(aes(fill=subtype),draw_quantiles = .5,)+
  stat_compare_means()+
  geom_jitter()+
  facet_wrap(~subtype)+
  scale_fill_manual(name = "Subtype",values=cluster_colors)+
  theme_classic()




