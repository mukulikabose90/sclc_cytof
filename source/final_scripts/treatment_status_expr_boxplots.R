################################################################################
# This script plots violin plots of the expression of all protein markers between
# naive CTCs and CTCs treated with SOC
################################################################################
source("source/sclc_cytof_functions.R")

set.seed(42)
################################################################################
# Read in data
################################################################################
ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

################################################################################
# Create plot dataframe
################################################################################
markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam","SLUG", "PD-L1", "p-YAP", "CD44", "CD24","E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist")

y <- assay(ctcs, "exprs")

df <- data.frame(t(y), colData(ctcs), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")

gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(ctcs)))

plot_df <- gg_df %>%
  dplyr::filter(antigen %in% markers_to_use)

plot_df$antigen <- factor(plot_df$antigen, levels = markers_to_use)

plot_df$treatment_status <- ifelse(plot_df$treatment_status == "naive", "Naive", "SOC")

plot_df$treatment_status <- factor(plot_df$treatment_status, levels = c("Naive", "SOC"))

# Remove post-tarla samples
plot_df <- plot_df %>%
  filter(tarla != "post" | is.na(tarla))

################################################################################
# Plot violin plots
################################################################################
stat.test <- plot_df %>%
  group_by(antigen) %>%
  wilcox_test(expression ~ treatment_status) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test


p <- ggviolin(plot_df, x="treatment_status" ,y="expression", fill="treatment_status", lwd=.3, outlier.size = .1,draw_quantiles =0.5)+
  facet_wrap(~antigen,nrow=3)+
  ylim(0,9)+
  labs(y="Expression",
       x= "")+
  scale_fill_manual(name = "Subtype",values=c("#E63946","#457B9D"))+
  theme(axis.title = element_text(size=20),
        axis.text.x = element_text(size=16,angle=45,hjust=1,vjust=1),
        strip.text = element_text(face = "bold", size=16), 
        strip.background = element_blank())+
  rremove("legend")


stat.test <- stat.test %>% add_xy_position(x = "treatment_status")
p <- p + stat_pvalue_manual(stat.test, label = "p.adj.signif",size=5)

p
################################################################################
# Save figure
################################################################################
tiff(glue("figures/treatment_status_expression_violinplot.tiff"), width=360,height=200, units = "mm", res=600)
print(p)
dev.off()
