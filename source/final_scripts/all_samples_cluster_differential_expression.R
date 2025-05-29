source("source/sclc_cytof_functions.R")

set.seed(42)
################################################################################
# Read in data
################################################################################

sce <- readRDS("data/cytof_objects/sclc_all_samples_with_clusters.rds")

################################################################################
# Read in CTC clusters
################################################################################

ctc_clusters <- readRDS("data/ctc_clusters.rds")

################################################################################
# Run differential expression 
################################################################################
metric_to_use <- "median"

df <- cytof_de(sce, method = "wilcox", metric = metric_to_use, ident = "new_clusters")

################################################################################
# Create DE barplots
################################################################################
plot_df <- df %>% 
  mutate(significance = factor(ifelse(p_adj < .05, "*", ""), levels=c("*",""))) %>% 
  mutate(protein = reorder_within(protein, logfc, ident_list)) %>% 
  mutate(star_x = ifelse(logfc > 0, logfc+.25,logfc-.5)) %>% 
  mutate(up_down= ifelse(logfc>0,"up","down"))

# Reorder cluster factor
plot_df$ident_list <- paste0("Cluster ",plot_df$ident_list)
plot_df$ident_list <- factor(plot_df$ident_list, levels=paste0("Cluster ", c(1:10)))

# Dataframe for plotting only CTC clusters
ctc_df <- plot_df %>% 
  dplyr::filter(ident_list %in% paste0("Cluster ", ctc_clusters))


x_axis_label <- gsub("m","M",metric_to_use)

# Create plot for CTC clusters
p1 <- ggplot(ctc_df,aes(x=as.numeric(logfc), y=protein, fill=logfc))+
  geom_col(color="darkgray",linewidth=.001)+
  geom_text(aes(x=star_x, label=significance), size = 3)+
  facet_wrap(~ident_list, scales = "free_y",nrow=1)+
  scale_y_reordered()+
  labs(fill = "")+
  guides(fill="none")+
  # scale_fill_manual(values=c("royalblue4", "firebrick"))+
  scale_fill_gradient2(low = "royalblue4", mid = "white", high = "firebrick", midpoint = 0)+
  ylab("Protein")+
  xlab(glue("{x_axis_label} log(FC)"))+
  theme_classic()+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=8), 
        axis.text = element_text(color = "black", size=6),
        axis.title = element_text(size=8),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8))

p1

# Create plot for all clusters
p2 <- ggplot(plot_df,aes(x=as.numeric(logfc), y=protein, fill=logfc))+
  geom_col(color="darkgray",size=.001)+
  geom_text(aes(x=star_x, label=significance), size = 3)+
  facet_wrap(~ident_list, scales = "free_y", nrow=3)+
  scale_y_reordered()+
  labs(fill = "")+
  guides(fill="none")+
  # scale_fill_manual(values=c("royalblue", "firebrick1"))+
  scale_fill_gradient2(low = "royalblue4", mid = "white", high = "firebrick", midpoint = 0)+
  ylab("Protein")+
  xlab(glue("{x_axis_label} log(FC)"))+
  theme_classic()+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=8), 
        axis.text = element_text(color = "black", size=6),
        axis.title = element_text(size=8),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8))

p2
################################################################################
# Save figures
################################################################################
jpeg(glue("figures/all_samples_cluster_{metric_to_use}_de_signif.jpg"), width=120,height=50, units = "mm", res=1000)
print(p1)
dev.off()

jpeg(glue("figures/all_samples_cluster_{metric_to_use}_de_all.jpg"), width=140,height=160, units = "mm", res=1000)
print(p2)
dev.off()
