################################################################################
# This script plots barplots displaying the proportion of cell of each subtype
# in every sample in the analysis. Split by treatment status
################################################################################
source("source/sclc_cytof_functions.R")

set.seed(42)
################################################################################
# Read in data
################################################################################
ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

################################################################################
# Set up plot dataframe
################################################################################
cluster_colors <- c("#dd4b33", "#F1FAEE", "#A8DADC", "#457B9D")

plot_df <-  ctcs@colData %>% 
  as.data.frame()  %>% 
  count(collection_id,subtype,treatment_status) %>% 
  group_by(collection_id) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n/total)*100) %>% 
  group_by(collection_id)

plot_df <- plot_df %>% 
  filter(total >= 10)

sample_order <- plot_df %>% 
  filter(subtype == "M") %>% 
  arrange(desc(freq)) %>% 
  pull(collection_id)

plot_df$collection_id <- factor(plot_df$collection_id, levels = sample_order)

plot_df$subtype <- factor(plot_df$subtype, levels = c("A","N","P","Mes"))

plot_df$treatment_status <- ifelse(plot_df$treatment_status == "naive", "Naive","Treated")
plot_df$treatment_status <- factor(plot_df$treatment_status, levels = c("Naive","Treated"))


p1 <- ggplot(plot_df)+
  geom_col(aes(x=collection_id,y=freq,fill=subtype))+
  geom_text(aes(label=total,x=collection_id,y=103),size=3)+
  # facet_wrap(~treatment_status,scales="free")+
  facet_grid(. ~ treatment_status, scales = "free", space = "free") +
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
################################################################################
# Save figure
################################################################################
tiff("figures/ctc_subtype_sample_proportion_barplots.tiff", width=300,height=100, units = "mm", res=600)
print(p1)
dev.off()

################################################################################
# Plot with tarla facet
################################################################################
plot_df <- ctcs@colData %>% 
  as.data.frame()

plot_df$treatment_status <- ifelse(plot_df$treatment_status == "naive", "Naive","CTX ± ICI")

plot_df$treatment_status <- ifelse(is.na(plot_df$tarla), plot_df$treatment_status,
                                   ifelse(plot_df$tarla == "pre", plot_df$treatment_status, "Tarla"))

plot_df <-  plot_df %>% 
  count(collection_id,subtype,treatment_status) %>% 
  group_by(collection_id) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n/total)*100) %>% 
  group_by(collection_id)


plot_df <- plot_df %>% 
  filter(total >= 10)

sample_order <- plot_df %>% 
  filter(subtype == "Mes") %>% 
  arrange(desc(freq)) %>% 
  pull(collection_id)

plot_df$collection_id <- factor(plot_df$collection_id, levels = sample_order)

plot_df$subtype <- factor(plot_df$subtype, levels = c("A","N","P","Mes"))

plot_df$treatment_status <- factor(plot_df$treatment_status, levels = c("Naive","CTX ± ICI","Tarla"))

plot_df <- plot_df %>% 
  filter(collection_id != "MDA-SC292-2")

plot_df$red_star <- ifelse(plot_df$total > 30, "*","")

p2 <- ggplot(plot_df)+
  geom_col(aes(x=collection_id,y=freq,fill=subtype))+
  geom_text(aes(label=total,x=collection_id,y=103),size=3)+
  geom_text(aes(label=red_star,x=collection_id,y=106),size=3,color="red")+
  # facet_wrap(~treatment_status,scales="free")+
  facet_grid(. ~ treatment_status, scales = "free", space = "free") +
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
################################################################################
# Save figure
################################################################################
tiff("figures/ctc_subtype_sample_proportion_barplots_tarla.tiff", width=300,height=100, units = "mm", res=600)
print(p2)
dev.off()

################################################################################
# Count
################################################################################

plot_df$patient_id <- sapply(as.character(plot_df$collection_id), function(x) strsplit(x,"-")[[1]][1])

plot_df %>% 
  pull(patient_id) %>% 
  unique() %>% 
  length()

plot_df %>% 
  pull(collection_id) %>% 
  unique() %>% 
  length()

# Naive patients
plot_df %>% 
  filter(treatment_status == "Naive") %>% 
  pull(patient_id) %>% 
  unique() %>% 
  length()

# Naive LBs
plot_df %>% 
  filter(treatment_status == "Naive") %>% 
  pull(collection_id) %>% 
  unique() %>% 
  length()

# SOC patients
plot_df %>% 
  filter(treatment_status == "SOC") %>% 
  pull(patient_id) %>% 
  unique() %>% 
  length()

# SOC LBs
plot_df %>% 
  filter(treatment_status == "SOC") %>% 
  pull(collection_id) %>% 
  unique() %>% 
  length()

# Post-tarla patients
plot_df %>% 
  filter(treatment_status == "Tarla") %>% 
  pull(patient_id) %>% 
  unique() %>% 
  length()

# Post-tarla LBs
plot_df %>% 
  filter(treatment_status == "Tarla") %>% 
  pull(collection_id) %>% 
  unique() %>% 
  length()
