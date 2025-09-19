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
  filter(total > 10)

sample_order <- plot_df %>% 
  filter(subtype == "I") %>% 
  arrange(desc(freq)) %>% 
  pull(collection_id)

plot_df$collection_id <- factor(plot_df$collection_id, levels = sample_order)

plot_df$subtype <- factor(plot_df$subtype, levels = c("A","N","P","I"))

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

################################################################################
# Save figure
################################################################################
tiff("figures/ctc_subtype_sample_proportion_barplots.tiff", width=300,height=100, units = "mm", res=600)
print(p1)
dev.off()
