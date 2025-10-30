################################################################################
# This script plots barplots displaying the proportion of cell of each subtype
# between pre-tarla CTCs and post-tarla CTCs 
################################################################################
source("source/sclc_cytof_functions.R")

set.seed(42)
################################################################################
# Read in data
################################################################################
ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

cluster_colors <- c("#dd4b33", "#F1FAEE", "#A8DADC", "#457B9D")

################################################################################
# Plot total subtype proportions between pre and post
################################################################################
curr_data <- ctcs

as.data.frame(curr_data@colData) %>% 
  filter(!is.na(tarla)) %>% 
  select(patient_id,tarla) %>% 
  distinct() %>% 
  count(patient_id) %>% 
  filter(n > 1)

plot_df <- as.data.frame(curr_data@colData) %>% 
  filter(!is.na(tarla)) %>% 
  select(tarla,subtype) %>% 
  dplyr::count(tarla,subtype) %>% 
  group_by(tarla) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n/total)*100) 

plot_df$total <- ifelse(plot_df$subtype == "Mes", plot_df$total,"")

plot_df$subtype <- factor(plot_df$subtype, levels=c("A","N","P","Mes"))

plot_df$tarla <- ifelse(plot_df$tarla == "pre","Pre-Tarlatamab","Post-Tarlatamab")
plot_df$tarla <- factor(plot_df$tarla, levels=c("Pre-Tarlatamab","Post-Tarlatamab"))

p <- ggplot(plot_df)+
  geom_col(aes(x=tarla,y=freq,fill=subtype))+
  geom_text(aes(label=total,x=tarla), y=105,size = 5)+
  ylim(0,105)+
  scale_fill_manual(name = "Subtype",values=cluster_colors)+
  labs(x="",
       y="Percentage",
       fill="Subtype")+
  theme_classic()+
  theme(axis.title = element_text(size=20),
        axis.text.x = element_text(size=12, color="black"),
        axis.text.y = element_text(size=16, color="black"),
        strip.text = element_text(face = "bold", size=12), 
        strip.background = element_blank(),
        legend.title = element_text(size=20),
        legend.text = element_text(size=18))


tiff("figures/tarla_subtype_barplots.tiff", width=140,height=100, units = "mm", res=600)
print(p)
dev.off()

################################################################################
# Plot subtype proportions for each patient (pre vs post)
################################################################################
# Select patients to use
patients_to_use <- ctcs %>%
  colData() %>%
  as.data.frame() %>%
  filter(!is.na(tarla)) %>%
  select(patient_id,tarla) %>%
  distinct() %>%
  count(patient_id) %>%
  filter(n > 1) %>%
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

patients_to_use <- patients_to_use[!patients_to_use %in% patients_to_remove]

curr_data <- ctcs[,ctcs$patient_id %in% patients_to_use]

plot_df <- as.data.frame(curr_data@colData) %>% 
  select(patient_id,tarla,subtype) %>% 
  dplyr::count(patient_id,tarla,subtype) %>% 
  group_by(patient_id,tarla) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n/total)*100) 

plot_df$patient_id <- as.character(plot_df$patient_id)

plot_df <- plot_df %>% 
  as.data.frame() %>% 
  tidyr::complete(patient_id,subtype,fill=list(n=0)) %>% 
  group_by(patient_id) %>% 
  mutate(count = max(total[which(!is.na(total))]))

plot_df <- plot_df %>% 
  dplyr::filter(!is.na(tarla))


plot_df$tarla <- ifelse(plot_df$tarla == "pre","Pre-Tarlatamab","Post-Tarlatamab")
plot_df$tarla <- factor(plot_df$tarla, levels=c("Pre-Tarlatamab","Post-Tarlatamab"))

plot_df$subtype <- factor(plot_df$subtype, levels=c("A","N","P","Mes"))

p <- ggplot(plot_df)+
  geom_col(aes(x=tarla,y=freq,fill=subtype))+
  geom_text(aes(label=total,x=tarla), y=105,size = 3)+
  facet_wrap(~patient_id, scales="free",ncol=5)+
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

tiff("figures/tarla_patient_subtype_barplots.tiff", width=250,height=100, units = "mm", res=1000)
print(p)
dev.off()

