source("source/sclc_cytof_functions.R")

set.seed(42)
################################################################################
# Read in data
################################################################################
ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

cluster_colors <- c("#dd4b33", "#F1FAEE", "#A8DADC", "#457B9D")
################################################################################
# Sample n cells from each patient
################################################################################
n_cells <- 30

sampled_data <- ctcs@colData %>%
  as.data.frame() %>%
  group_by(collection_id) %>%
  filter(n() >= n_cells) %>%        # keep only patients with ≥ n_cells cells
  slice_sample(n = n_cells) %>%
  ungroup()

lbs_to_use <- sampled_data$collection_id %>% 
  unique()

sampled_data$treatment_status <- ifelse(sampled_data$treatment_status == "naive","Naive","SOC")

sampled_data$treatment_status <- ifelse(is.na(sampled_data$tarla), sampled_data$treatment_status,
                                        ifelse(sampled_data$tarla == "pre", sampled_data$treatment_status, "Tarla"))


plot_df <-  sampled_data %>% 
  count(collection_id,subtype,treatment_status) %>% 
  group_by(collection_id) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n/total)*100) %>% 
  group_by(collection_id)


sample_order <- plot_df %>% 
  filter(subtype == "I") %>% 
  arrange(desc(freq)) %>% 
  pull(collection_id)


plot_df$subtype <- factor(plot_df$subtype, levels = c("A","N","P","Mes"))


p1 <- ggplot(plot_df)+
  geom_col(aes(x=collection_id,y=freq,fill=subtype))+
  geom_text(aes(label=total,x=collection_id,y=103),size=3)+
  # facet_wrap(~treatment_status,scales="free")+
  facet_grid(. ~ treatment_status, scales = "free", space = "free",labeller = as_labeller(c("Naive"="Naive",
                                                                                            "SOC"="CTX ± ICI",
                                                                                            "Tarla"="Tarla"))) +
  scale_fill_manual(values = cluster_colors)+
  theme_classic()+
  labs(x="Sample",
       y="Percentage",
       fill="Subtype")+
  theme(axis.text.x = element_text(size=12, angle=45,hjust=1,vjust=1),
        axis.text.y = element_text(size=16),
        axis.title = element_text(size=20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=18),
        strip.text = element_text(size=20))

p1

tiff("figures/downsampled_ctc_subtype_sample_proportion_barplots_tarla.tiff", width=300,height=100, units = "mm", res=600)
print(p1)
dev.off()

################################################################################

sampled_data <- ctcs@colData %>% 
  as.data.frame() %>% 
  filter(collection_id %in% lbs_to_use)

sampled_data$treatment_status <- ifelse(sampled_data$treatment_status == "naive","Naive","SOC")

sampled_data$treatment_status <- ifelse(is.na(sampled_data$tarla), sampled_data$treatment_status,
                                        ifelse(sampled_data$tarla == "pre", sampled_data$treatment_status, "Tarla"))


plot_df <-  sampled_data %>% 
  count(collection_id,subtype,treatment_status) %>% 
  group_by(collection_id) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n/total)*100) %>% 
  group_by(collection_id)


sample_order <- plot_df %>% 
  filter(subtype == "I") %>% 
  arrange(desc(freq)) %>% 
  pull(collection_id)

plot_df$subtype <- factor(plot_df$subtype, levels = c("A","N","P","Mes"))


p2 <- ggplot(plot_df)+
  geom_col(aes(x=collection_id,y=freq,fill=subtype))+
  geom_text(aes(label=total,x=collection_id,y=103),size=3)+
  # facet_wrap(~treatment_status,scales="free")+
  facet_grid(. ~ treatment_status, scales = "free", space = "free",labeller = as_labeller(c("Naive"="Naive",
                                                                                            "SOC"="CTX ± ICI",
                                                                                            "Tarla"="Tarla"))) +
  scale_fill_manual(values = cluster_colors)+
  theme_classic()+
  labs(x="Sample",
       y="Percentage",
       fill="Subtype")+
  theme(axis.text.x = element_text(size=12, angle=45,hjust=1,vjust=1),
        axis.text.y = element_text(size=16),
        axis.title = element_text(size=20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=18),
        strip.text = element_text(size=20))

p2

tiff("figures/downsampled_all_cells_ctc_subtype_sample_proportion_barplots_tarla.tiff", width=300,height=100, units = "mm", res=600)
print(p2)
dev.off()

################################################################################

sampled_data <- ctcs@colData %>% 
  as.data.frame() %>% 
  filter(collection_id %in% lbs_to_use)




sampled_data$treatment_status <- ifelse(sampled_data$treatment_status == "naive","Naive","SOC")

sampled_data$treatment_status <- ifelse(is.na(sampled_data$tarla), sampled_data$treatment_status,
                                        ifelse(sampled_data$tarla == "pre", sampled_data$treatment_status, "Tarla"))


matched_patients <- sampled_data %>% 
  count(patient_id, treatment_status) %>% 
  count(patient_id) %>% 
  filter(n > 1) %>% 
  pull(patient_id)


plot_df <-  sampled_data %>% 
  filter(patient_id %in% matched_patients) %>% 
  count(patient_id,subtype,treatment_status) %>% 
  group_by(patient_id,treatment_status) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n/total)*100) %>% 
  group_by(patient_id)


plot_df$subtype <- factor(plot_df$subtype, levels = c("A","N","P","Mes"))

p3 <- ggplot(plot_df)+
  geom_col(aes(x=treatment_status,y=freq,fill=subtype))+
  geom_text(aes(label=total,x=treatment_status,y=103),size=3)+
  # facet_wrap(~treatment_status,scales="free")+
  facet_grid(. ~ patient_id, scales = "free", space = "free") +
  scale_fill_manual(values = cluster_colors)+
  theme_classic()+
  labs(x="Sample",
       y="Percentage",
       fill="Subtype")+
  theme(axis.text.x = element_text(size=12, angle=45,hjust=1,vjust=1),
        axis.text.y = element_text(size=16),
        axis.title = element_text(size=20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=18),
        strip.text = element_text(size=20))

p3
