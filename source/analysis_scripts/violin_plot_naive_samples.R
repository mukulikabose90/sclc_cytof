source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam", "E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist", "SLUG", "PD-L1", "p-YAP", "CD44", "CD24")
################################################################################

################################################################################
sce <- ctcs

y <- assay(sce, "exprs")

df <- data.frame(t(y), colData(sce), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")

gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(sce)))

plot_df <- gg_df %>%
  dplyr::filter(antigen %in% markers_to_use)

plot_df$antigen <- factor(plot_df$antigen, levels = markers_to_use)

# plot_df$treatment_status <- ifelse(plot_df$treatment_status == "naive", "Naive", "Treated")
# 
# plot_df$treatment_status <- factor(plot_df$treatment_status, levels = c("Naive", "Treated"))

################################################################################

# plot_df <- plot_df %>% 
#   filter(treatment_status =="naive" & collection_id %in% c("SC506-1","SC293-1"))

plot_df <- plot_df %>% 
  filter(treatment_status =="naive")


p <- ggviolin(plot_df, x="collection_id" ,y="expression", fill="collection_id", lwd=.3, outlier.size = .1,draw_quantiles =0.5)+
  stat_compare_means(aes(group = treatment_status), label = "p.signif", label.x.npc = "center", label.y = 5.5,size=4.5)+ 
  facet_wrap(~antigen,nrow=3)+
  ylim(0,NA)+
  labs(y="Expression",
       x= "")+
  # scale_fill_manual(name = "Sample",values=c("#E63946","#457B9D"))+
  theme(axis.title = element_text(size=14),
        axis.text.x = element_text(size=10,angle=45,hjust=1,vjust=1),
        strip.text = element_text(face = "bold", size=12), 
        strip.background = element_blank())+
  rremove("legend")



p
