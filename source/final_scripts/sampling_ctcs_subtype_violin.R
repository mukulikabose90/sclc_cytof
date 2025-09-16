

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

patients_to_use <- as.data.frame(curr_data@colData) %>% 
  filter(tarla != "post" | is.na(tarla)) %>% 
  count(collection_id) %>% 
  filter(n >= 10) %>% 
  pull(collection_id) %>% 
  as.character()

as.data.frame(curr_data@colData) %>% 
  filter(tarla != "post" | is.na(tarla)) %>% 
  filter(patient_id %in% patients_to_use)

sampled_results <- list()
for(i in 1:100){
  results_list <- list()
  for(curr_patient in patients_to_use) {
    
    results_list[[curr_patient]] <- as.data.frame(curr_data@colData) %>% 
      filter(collection_id == curr_patient) %>% 
      sample_n(., 10)
    
    
    
    
  }
  
  sampled_df <- bind_rows(results_list)
  
  plot_df <- sampled_df %>% 
    dplyr::count(treatment_status,subtype) %>% 
    group_by(treatment_status) %>% 
    mutate(total = sum(n)) %>% 
    mutate(freq = (n/total)*100) 
  
  sampled_results[[i]] <- plot_df
}

final_df <- bind_rows(sampled_results)

ggplot(final_df,aes(x=treatment_status,y=freq))+
  geom_violin(aes(fill=subtype),draw_quantiles = .5)+
  stat_compare_means()+
  facet_wrap(~subtype)




plot_df$subtype <- factor(plot_df$subtype, levels=c("A","N","P",'I'))

plot_df$treatment_status <- ifelse(plot_df$treatment_status == "naive","Naive","Treated")


p <- ggplot(plot_df)+
  geom_col(aes(x=treatment_status,y=freq,fill=subtype))+
  geom_text(aes(label=total,x=treatment_status), y=105,size = 5)+
  ylim(0,105)+
  scale_fill_manual(name = "Subtype",values=cluster_colors)+
  labs(x="",
       y="Percentage",
       fill="Subtype")+
  theme_classic()+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=16, color="black"),
        strip.text = element_text(face = "bold", size=12), 
        strip.background = element_blank(),
        legend.title = element_text(size=20),
        legend.text = element_text(size=18))

p





