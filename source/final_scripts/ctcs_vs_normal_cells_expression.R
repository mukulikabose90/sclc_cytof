################################################################################
# This script plots the expression of p-Rb in:
#   - Normal cells vs CTCs
#   - CTCs vs non-CTCS (from cancer enriched clusters)
################################################################################
source("source/sclc_cytof_functions.R")

set.seed(42)
################################################################################
# Read in data
################################################################################
# CTCs data
ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

# All cells data
sce <- readRDS("data/cytof_objects/sclc_all_samples_object.rds")

sce$ctc <- ifelse(sce$cell_id %in% ctcs$cell_id, "ctc","normal")

# Subset to only cells from normal samples and CTCs
sce <- sce[,sce$condition == "normal" | sce$ctc == "ctc"]

################################################################################
# Plot CTCs vs normal cells
################################################################################

markers_to_use <- c("p-Rb")

y <- assay(sce, "exprs")

df <- data.frame(t(y), colData(sce), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")

gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
               id.vars = names(colData(sce)))


plot_df <- gg_df %>%
  dplyr::filter(antigen %in% markers_to_use)

plot_df$ctc <- ifelse(plot_df$ctc == "normal", "Normal","CTCs")
plot_df$ctc <- factor(plot_df$ctc, levels=c("CTCs","Normal"))

p1 <- ggviolin(plot_df, x="ctc" ,y="expression", fill="ctc",draw_quantiles = 0.5)+
  stat_compare_means(comparisons = list(c("CTCs","Normal")),label.y = 6,tip.length = 0)+
  facet_wrap(~antigen)+
  ylim(0,6.5)+
  labs(y="Expression",
       x= "")+
  scale_fill_manual(name = "Subtype",values=c("#E63946","#457B9D"))+
  theme(axis.title = element_text(size=14),
        axis.text.x = element_text(size=10,angle=45,hjust=1,vjust=1),
        strip.text = element_text(face = "bold", size=12), 
        strip.background = element_blank())+
  rremove("legend")
p1
################################################################################
# Plot CTCs vs non-CTCs in cancer enriched clusters
################################################################################
# Read in cancer enriched clusters data
sce <- readRDS("data/cytof_objects/cancer_enriched_with_clusters.rds")

sce$ctc <- ifelse(sce$new_clusters %in% c(2,4,5,6,7,8), "CTCs","Non-CTCs")

markers_to_use <- c("p-Rb")

y <- assay(sce, "exprs")

df <- data.frame(t(y), colData(sce), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")

gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(sce)))


plot_df <- gg_df %>%
  dplyr::filter(antigen %in% markers_to_use)

plot_df$ctc <- factor(plot_df$ctc, levels=c("CTCs","Non-CTCs"))


p2 <- ggviolin(plot_df, x="ctc" ,y="expression", fill="ctc",draw_quantiles = 0.5)+
  stat_compare_means(comparisons = list(c("CTCs","Non-CTCs")), label.y = 5.2,tip.length = 0)+
  facet_wrap(~antigen)+
  ylim(0,6.5)+
  labs(y="Expression",
       x= "")+
  scale_fill_manual(name = "Subtype",values=c("#E63946","#457B9D"))+
  theme(axis.title = element_text(size=14),
        axis.text.x = element_text(size=10,angle=45,hjust=1,vjust=1),
        strip.text = element_text(face = "bold", size=12), 
        strip.background = element_blank())+
  rremove("legend")

################################################################################
# Save figures
################################################################################

tiff(glue("figures/ctcs_vs_normal_rb_violin_plot.tiff"), width=100,height=100, units = "mm", res=1000)
print(p1)
dev.off()

tiff(glue("figures/ctcs_vs_nonctcs_rb_violin_plot.tiff"), width=100,height=100, units = "mm", res=1000)
print(p2)
dev.off()

