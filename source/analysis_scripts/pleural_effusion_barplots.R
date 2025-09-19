source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

cluster_colors <- c("#dd4b33", "#F1FAEE", "#A8DADC", "#457B9D")

################################################################################
# Plot total subtype proportions between pre and post
################################################################################
curr_data <- ctcs

samples_to_use <- c("SC506-1P","SC443-1P","SC399-3P")
samples_to_use <- c("SC443")


plot_df <- as.data.frame(curr_data@colData) %>% 
  filter(patient_id == samples_to_use) %>% 
  select(collection_id,subtype) %>% 
  dplyr::count(collection_id,subtype) %>% 
  group_by(collection_id) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n/total)*100) 

plot_df$total <- ifelse(plot_df$subtype == "I", plot_df$total,"")

plot_df$subtype <- factor(plot_df$subtype, levels=c("A","N","P",'I'))





p <- ggplot(plot_df)+
  geom_col(aes(x=collection_id,y=freq,fill=subtype))+
  geom_text(aes(label=total,x=collection_id), y=105,size = 3)+
  ylim(0,105)+
  scale_fill_manual(name = "Subtype",values=cluster_colors)+
  labs(x="",
       y="Percentage",
       fill="Subtype")+
  theme_classic()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=10, color="black"),
        strip.text = element_text(face = "bold", size=12), 
        strip.background = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12))

p

jpeg(glue("figures/pleural_effusion_subtype_barplots.jpg"), width=260,height=150, units = "mm", res=1000)
print(p)
dev.off()
