################################################################################
# This script plots violin plots for protein expression in patients that had PDXs
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
markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam", "E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist", "SLUG", "PD-L1", "p-YAP", "CD44", "CD24")

y <- assay(ctcs, "exprs")

df <- data.frame(t(y), colData(ctcs), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")

gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(ctcs)))

plot_df <- gg_df %>%
  dplyr::filter(antigen %in% markers_to_use)

plot_df$antigen <- factor(plot_df$antigen, levels = markers_to_use)

plot_df$treatment_status <- ifelse(plot_df$treatment_status == "naive", "Naive", "Treated")

plot_df$treatment_status <- factor(plot_df$treatment_status, levels = c("Naive", "Treated"))

patients_to_use <- c("SC293","SC506","SC443")

plot_df <- plot_df %>%
  filter(patient_id %in% patients_to_use & treatment_status == "Naive")

################################################################################
# Plot violin plots
################################################################################

my_comparisons <- list(c("SC293", "SC506"), c("SC506", "SC443"), c("SC293", "SC443"))

stat.test <- plot_df %>%
  group_by(antigen) %>%
  wilcox_test(expression ~ patient_id, comparisons = my_comparisons) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

plot_df$patient_id <- factor(plot_df$patient_id, levels=c("SC293","SC506","SC443"))

p <- ggviolin(plot_df, x="patient_id" ,y="expression", fill="patient_id", lwd=.3, outlier.size = .1,draw_quantiles =0.5)+
  facet_wrap(~antigen,nrow=2)+
  ylim(0,NA)+
  labs(y="Expression",
       x= "")+
  scale_fill_manual(name = "Patient",values=c("#E63946","#457B9D","#956789"))+
  theme(axis.title = element_text(size=14),
        axis.text.x = element_text(size=10,angle=45,hjust=1,vjust=1),
        strip.text = element_text(face = "bold", size=12), 
        strip.background = element_blank())+
  rremove("legend")

stat.test <- stat.test %>% 
  add_xy_position(x = "patient_id")

p <- p + stat_pvalue_manual(stat.test, label = "p.adj.signif",xmin = "group1")

################################################################################
# Save figure
################################################################################
tiff(glue("figures/xenograft_samples_expression_violinplot.tiff"), width=360,height=140, units = "mm", res=1000)
print(p)
dev.off()
