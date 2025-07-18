source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

cluster_colors <- c("#dd4b33", "#F1FAEE", "#A8DADC", "#457B9D")

################################################################################
# Plot subtype proportions for each patient across all samples

long_patients <- as.data.frame(ctcs@colData) %>% 
  filter(tarla != "post" | is.na(tarla)) %>% 
  select(patient_id,sample_num) %>% 
  distinct() %>% 
  dplyr::count(patient_id) %>% 
  dplyr::filter(n > 1) %>% 
  pull(patient_id) %>% 
  as.character()

# Remove samples with too few cells
patients_to_remove <- ctcs %>% 
  colData() %>% 
  as.data.frame() %>% 
  count(patient_id,sample_num) %>% 
  filter(n < 10) %>% 
  pull(patient_id) %>% 
  unique() %>% 
  as.character()

long_patients <- long_patients[!long_patients %in% patients_to_remove]

long_data <- ctcs[,ctcs$patient_id %in% long_patients]

#Calculate proportions of each subtype for each sample
plot_df <- as.data.frame(long_data@colData) %>% 
  select(patient_id,sample_num,subtype) %>% 
  dplyr::count(patient_id,sample_num,subtype) %>% 
  group_by(patient_id,sample_num) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n/total)*100) %>% 
  mutate(sample_id = paste0(patient_id,"-",sample_num)) %>% 
  as.data.frame() %>% 
  tidyr::complete(sample_id,subtype,fill=list(n=0)) %>% 
  group_by(sample_id) %>% 
  mutate(count = max(total[which(!is.na(total))])) %>% 
  dplyr::filter(!is.na(patient_id))

#Rearrange subtypes
plot_df$subtype <- factor(plot_df$subtype, levels=c("A","N","P",'I'))

#Change sample number to sample order
df_renumbered <- plot_df %>%
  select(patient_id,sample_num) %>% 
  distinct() %>% 
  group_by(patient_id) %>%
  arrange(sample_num, .by_group = TRUE) %>%
  mutate(new_sample_id = row_number()) %>%
  ungroup() %>% 
  select(sample_id,new_sample_id)


plot_df <- merge(plot_df,df_renumbered,by="sample_id")

plot_df$new_sample_id <- factor(plot_df$new_sample_id)

p <- ggplot(plot_df)+
  geom_col(aes(x=new_sample_id,y=freq,fill=subtype))+
  geom_text(aes(label=total,x=new_sample_id), y=105,size = 3)+
  facet_wrap(~patient_id, scales="free",ncol=5)+
  ylim(0,105)+
  scale_fill_manual(name = "Subtype",values=cluster_colors)+
  labs(x="Sample Number",
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


jpeg("figures/all_patients_longitudinal_subtype_barplots.jpg", width=225,height=150, units = "mm", res=1000)
print(p)
dev.off()








