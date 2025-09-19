################################################################################
# This script plots the subtype proportion barplots for the patients that had
# PDXs 
################################################################################
source("source/sclc_cytof_functions.R")

set.seed(42)
################################################################################
# Read in data
################################################################################
ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

################################################################################
# Plot total subtype proportions between pre and post
################################################################################
curr_data <- ctcs

cluster_colors <- c("#dd4b33", "#F1FAEE", "#A8DADC", "#457B9D")

patients_to_use <- c("SC293","SC506","SC443")

plot_df <- as.data.frame(curr_data@colData) %>% 
  filter(patient_id %in% patients_to_use & treatment_status == "naive") %>% 
  select(patient_id,subtype) %>% 
  dplyr::count(patient_id,subtype) %>% 
  group_by(patient_id) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n/total)*100) 

plot_df$total <- ifelse(plot_df$subtype == "I", plot_df$total,"")

plot_df$subtype <- factor(plot_df$subtype, levels=c("A","N","P",'I'))

plot_df$patient_id <- factor(plot_df$patient_id, levels=c("SC293","SC506","SC443"))

################################################################################
# Plot barplots
################################################################################

p <- ggplot(plot_df)+
  geom_col(aes(x=patient_id,y=freq,fill=subtype))+
  geom_text(aes(label=total,x=patient_id), y=105,size = 5)+
  ylim(0,105)+
  scale_fill_manual(name = "Subtype",values=cluster_colors)+
  labs(x="",
       y="Percentage",
       fill="Subtype")+
  theme_classic()+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=12, color="black"),
        strip.text = element_text(face = "bold", size=12), 
        strip.background = element_blank(),
        legend.title = element_text(size=20),
        legend.text = element_text(size=18))

################################################################################
# Save figure
################################################################################
tiff(glue("figures/xenograft_subtype_barplots.tiff"), width=140,height=100, units = "mm", res=600)
print(p)
dev.off()
