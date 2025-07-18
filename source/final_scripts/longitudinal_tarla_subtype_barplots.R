source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

################################################################################
# Plot subtype proportions for each patient (pre vs post tarla)

long_patients <- as.data.frame(ctcs@colData) %>% 
  filter(!is.na(tarla)) %>%
  select(patient_id,tarla) %>% 
  distinct() %>% 
  dplyr::count(patient_id) %>% 
  dplyr::filter(n > 1) %>% 
  pull(patient_id) %>% 
  as.character()

long_data <- ctcs[,ctcs$patient_id %in% long_patients]

plot_df <- as.data.frame(long_data@colData) %>% 
  select(patient_id,tarla,subtype) %>% 
  dplyr::count(patient_id,tarla,subtype) %>% 
  group_by(patient_id,tarla) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = n/total) 

plot_df$patient_id <- as.character(plot_df$patient_id)

plot_df <- plot_df %>% 
  dplyr::filter(!is.na(patient_id))

plot_df$freq <- plot_df$freq*100

plot_df$subtype <- factor(plot_df$subtype, levels=c("A","N","P",'I'))

plot_df$tarla <- ifelse(plot_df$tarla == "pre","Pre-Tarlatamab","Post-Tarlatamab")

plot_df$tarla <- factor(plot_df$tarla, levels=c("Pre-Tarlatamab","Post-Tarlatamab"))

p <- ggplot(plot_df)+
  geom_col(aes(x=tarla,y=freq,fill=subtype))+
  geom_text(aes(label=total,x=tarla), y=105,size = 3)+
  facet_wrap(~patient_id, scales="free",ncol=2)+
  ylim(0,105)+
  scale_fill_manual(name = "Subtype",values=cluster_colors)+
  labs(x="",
       y="Percentage",
       fill="Subtype")+
  theme_classic()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=8, color="black"),
        strip.text = element_text(face = "bold", size=12), 
        strip.background = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12))


p

jpeg("figures/longitudinal_tarla_barplots.jpg", width=150,height=100, units = "mm", res=1000)
print(p)
dev.off()

################################################################################
# Plot total subtype proportions between pre and post tarla

plot_df <- as.data.frame(ctcs@colData) %>% 
  filter(!is.na(tarla)) %>% 
  select(tarla,subtype) %>% 
  dplyr::count(tarla,subtype) %>% 
  group_by(tarla) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = n/total) 

plot_df$total <- ifelse(plot_df$subtype == "I", plot_df$total,"")

plot_df$freq <- plot_df$freq*100

plot_df$tarla <- ifelse(plot_df$tarla == "pre","Pre-Tarlatamab","Post-Tarlatamab")

plot_df$tarla <- factor(plot_df$tarla, levels=c("Pre-Tarlatamab","Post-Tarlatamab"))

plot_df$subtype <- factor(plot_df$subtype, levels=c("A","N","P",'I'))

p <- ggplot(plot_df)+
  geom_col(aes(x=tarla,y=freq,fill=subtype))+
  geom_text(aes(label=total,x=tarla), y=105,size = 3)+
  ylim(0,105)+
  scale_fill_manual(name = "Subtype",values=cluster_colors)+
  labs(x="",
       y="Percentage",
       fill="Subtype")+
  theme_classic()+
  theme(axis.title = element_text(size=14,color="black"),
        axis.text = element_text(size=10, color="black"),
        strip.text = element_text(face = "bold", size=12), 
        strip.background = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12))


p


jpeg("figures/tarla_subtype_barplots.jpg", width=120,height=100, units = "mm", res=1000)
print(p)
dev.off()

################################################################################
