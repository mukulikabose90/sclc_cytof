source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################

markers_to_use <- c("p-YAP")

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

curr_data <- ctcs[,ctcs$patient_id %in% c("SC293","SC338","SC454")]

y <- assay(curr_data, "exprs")

df <- data.frame(t(y), colData(curr_data), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")

gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(curr_data)))


plot_df <- gg_df %>%
  dplyr::filter(antigen %in% markers_to_use)

plot_df$antigen <- factor(plot_df$antigen, levels = markers_to_use)


plot_df$sample_num <- ifelse(plot_df$sample_num == "4" & plot_df$patient_id == "SC338", "1", plot_df$sample_num)
plot_df$sample_num <- ifelse(plot_df$sample_num == "1", "1", "2")

################################################################################

stat.test <- plot_df %>%
  group_by(patient_id) %>%
  wilcox_test(expression ~ sample_num) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

plot_df$patient_id <- factor(plot_df$patient_id,levels = c("SC293","SC338","SC454"))

p <- ggviolin(plot_df, x="patient_id" ,y="expression", fill="sample_num", lwd=.3,
              outlier.size = .1, draw_quantiles = .5, position = position_dodge(width = 1))+
  facet_wrap(~antigen,nrow=2)+
  ylim(0,NA)+
  labs(y="Expression",
       x= "")+
  scale_fill_manual(name = "Sample Number",values=c("#E63946","#457B9D"))+
  theme(axis.title = element_text(size=20),
        axis.text.x = element_text(size=16,angle=45,hjust=1,vjust=1),
        strip.text = element_text(face = "bold", size=16), 
        strip.background = element_blank(),
        legend.position = "right")
  rremove("legend")

p

stat.test <- stat.test %>% add_xy_position(x="patient_id")

p <- p + stat_pvalue_manual(stat.test, label = "p.adj.signif", size=5, tip.length = 0,
                            y.position = c(4,2,2.2))

tiff(glue("figures/pyap_longitudinal_expression_violin.tiff"), width=260,height=150, units = "mm", res=600)
print(p)
dev.off()

plot_df %>% 
  count(patient_id,sample_num)







