source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################
sce <- readRDS("data/cytof_objects/sclc_all_samples_with_clusters.rds")

ctc_clusters <- readRDS("data/ctc_clusters.rds")

# sce$ctc <- ifelse(sce$new_clusters %in% ctc_clusters & sce$condition == "cancer", "ctc", "non-ctc")
sce$ctc <- ifelse(sce$new_clusters %in% ctc_clusters, "ctc", "non-ctc")

y <- assay(sce, "exprs")

df <- data.frame(t(y), colData(sce), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")

gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(sce)))

sclc_tfs <- c("NeuroD1","ASCL1","POU2F3","p-Rb")

# sclc_tfs <- markers

temp <- gg_df %>%
  dplyr::filter(antigen %in% sclc_tfs)

temp$ctc <- factor(temp$ctc, levels=c("ctc","non-ctc"))


p <- ggboxplot(temp, x="new_clusters",y="expression", fill="ctc")

p1 <- p + facet_wrap(~antigen, scales = "free_y")+
  stat_compare_means()+
  scale_fill_manual(values=c("#E57373","#64B5F6"))

p <- ggboxplot(temp, x="ctc",y="expression", fill="ctc")

p2 <- p + facet_wrap(~antigen, scales = "free_y")+
  stat_compare_means()+
  scale_fill_manual(values=c("#E57373","#64B5F6"))



jpeg(glue("figures/all_samples_cluster_tf_boxplot.jpg"), width=160,height=160, units = "mm", res=1000)
print(p1)
dev.off()

jpeg(glue("figures/all_samples_ctc_tf_boxplot.jpg"), width=160,height=160, units = "mm", res=1000)
print(p2)
dev.off()
