source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam", "E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist", "SLUG", "PD-L1", "p-YAP", "CD44", "CD24")

################################################################################

sce <- ctcs

y <- assay(sce, "exprs")

df <- data.frame(t(y), colData(sce), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")

gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(sce)))


plot_df <- gg_df %>%
  dplyr::filter(antigen %in% markers_to_use)

unique(plot_df$tarla)

# Remove non-tarla samples
plot_df <- plot_df %>%
  filter(!is.na(tarla))

plot_df$antigen <- factor(plot_df$antigen, levels = markers_to_use)

plot_df$tarla <- ifelse(plot_df$tarla == "pre", "Pre", "Post")

plot_df$tarla <- factor(plot_df$tarla, levels=c("Pre", "Post"))

################################################################################

stat.test <- plot_df %>%
  group_by(antigen) %>%
  wilcox_test(expression ~ tarla) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

p <- ggviolin(plot_df, x="tarla" ,y="expression", fill="tarla", lwd=.3, outlier.size = .1, draw_quantiles = .5)+
  # stat_compare_means(aes(group = tarla), label = "p.signif", label.x.npc = "center", label.y = 5.5,size=4.5)+ 
  facet_wrap(~antigen,nrow=2)+
  ylim(0,9)+
  labs(y="Expression",
       x= "Tarlatamab Status")+
  scale_fill_manual(name = "Subtype",values=c("#E63946","#457B9D"))+
  theme(axis.title = element_text(size=14),
        axis.text.x = element_text(size=10,angle=45,hjust=1,vjust=1),
        strip.text = element_text(face = "bold", size=12), 
        strip.background = element_blank())+
  rremove("legend")

stat.test <- stat.test %>% add_xy_position(x = "tarla")
p <- p + stat_pvalue_manual(stat.test, label = "p.adj.signif")


jpeg(glue("figures/pre_post_tarla_expression_violinplot.jpg"), width=360,height=140, units = "mm", res=1000)
print(p)
dev.off()

