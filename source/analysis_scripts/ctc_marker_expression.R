source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################
library(ggbeeswarm)
ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

y <- assay(ctcs, "exprs")

df <- data.frame(t(y), colData(ctcs), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")

gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(ctcs)))


unique(gg_df$antigen)

markers_to_use <- c("PD-L1","DLL3","MUC-1","p-YAP","Vimentin","E-Cad","EpCAM")
markers_to_use <- c("Vimentin","E-Cad","EpCAM")

markers_to_use <- c("PD-L1","DLL3")
markers_to_use <- c("ASCL1","Twist","SLUG")


markers_to_use <- c("Ir191","Ir193")


plot_df <- gg_df %>%
  dplyr::filter(antigen %in% markers_to_use)



ggboxplot(plot_df, x="ctc",y="expression",fill="ctc",add = "jitter", add.params = list(size = 1, alpha = 0.025))+
  facet_wrap(~antigen,scales="free")+
  stat_compare_means(label = "p.signif",  label.y=3,tip.length = 0, comparisons = list(c("naive","treated")))+
  labs(y=glue("Expression"),
       x="")+
  rremove("legend")+
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=20))


plot_df <- plot_df %>% 
  filter(!is.na(tarla))


plot_df$tarla <- factor(plot_df$tarla, levels=c("pre","post"))

ggboxplot(plot_df, x="tarla",y="expression",fill="tarla",add = "jitter", add.params = list(size = 1, alpha = 0.025))+
  facet_wrap(~antigen)+
  stat_compare_means(label = "p.signif",  label.y=2, tip.length = 0, size=8, comparisons = list(c("pre","post")))+
  labs(y=glue("Expression"),
       x="")+
  ylim(c(0,2.5))+
  rremove("legend")+
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=20))
