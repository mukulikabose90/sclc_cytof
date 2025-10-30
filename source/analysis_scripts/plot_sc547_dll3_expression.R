
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


curr_data <- ctcs[,ctcs$patient_id == "MDA-SC547"]
################################################################################
# Create plot dataframe
################################################################################
# markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam","SLUG", "PD-L1", "p-YAP", "CD44", "CD24","E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist")
markers_to_use <- c("DLL3")

y <- assay(curr_data, "exprs")

df <- data.frame(t(y), colData(curr_data), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")

gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(curr_data)))

plot_df <- gg_df %>%
  dplyr::filter(antigen %in% markers_to_use)

plot_df$antigen <- factor(plot_df$antigen, levels = markers_to_use)



plot_df$treatment_status <- ifelse(is.na(plot_df$tarla), plot_df$treatment_status,
                                   ifelse(plot_df$tarla == "pre", plot_df$treatment_status, "Tarla"))

# plot_df$treatment_status <- ifelse(plot_df$treatment_status == "naive", "Naive", "SOC")

plot_df$treatment_status <- factor(plot_df$treatment_status, levels = c("Naive", "CTX + ICI","Tarla"))

# Remove post-tarla samples
# plot_df <- plot_df %>%
#   filter(tarla != "post" | is.na(tarla))

################################################################################
# Plot violin plots
################################################################################

plot_df$collection_id
stat.test <- plot_df %>%
  wilcox_test(expression ~ collection_id, comparisons = list(c("MDA-SC547-1","MDA-SC547-2"),c("MDA-SC547-1","MDA-SC547-3"),c("MDA-SC547-2","MDA-SC547-3"))) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance() %>% 
  add_xy_position(x = "collection_id")

stat.test

p1 <- ggviolin(plot_df, x="collection_id" ,y="expression", fill="collection_id", lwd=.3, outlier.size = .1,draw_quantiles =0.5,)+
  stat_pvalue_manual(stat.test, y.position = c(6.2,7,6.5),label = "p.adj.signif",size=5,tip.length = 0)+
  facet_wrap(~antigen,nrow=2)+
  ylim(0,NA)+
  labs(y="Expression",
       x= "")+
  scale_fill_manual(name = "Treatment",values=c("#1f77b3","#fe7e08","#2ba02d"))+
  theme(axis.title = element_text(size=20),
        axis.text.x = element_text(size=16,angle=45,hjust=1,vjust=1),
        axis.text.y = element_text(size=16),
        strip.text = element_text(face = "bold", size=16), 
        strip.background = element_blank(),
        legend.position = "none")


tiff(glue("figures/SC547_DLL3_expression_violinplot.tiff"), width=150,height=165, units = "mm", res=600)
print(p1)
dev.off()


################################################################################
# Plot SC547 with normals
################################################################################
ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")
ctc_547_data <- ctcs[,ctcs$patient_id == "MDA-SC547"]
curr_ctcs_ids <- ctc_547_data$cell_id

sce <- readRDS("data/cytof_objects/sclc_all_samples_object.rds")

curr_data <- sce[,sce$patient_id %in% paste0("NORMAL", 21:24) | sce$cell_id %in% curr_ctcs_ids]

dim(curr_data)

curr_data$collection_id <- paste0("MDA-",curr_data$collection_id)

################################################################################
# Subset to cancer cells in cancer_enriched cluster
################################################################################

markers_to_use <- c("DLL3")

y <- assay(curr_data, "exprs")

df <- data.frame(t(y), colData(curr_data), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")

gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(curr_data)))

plot_df <- gg_df %>%
  dplyr::filter(antigen %in% markers_to_use)

plot_df$antigen <- factor(plot_df$antigen, levels = markers_to_use)



plot_df$treatment_status <- ifelse(is.na(plot_df$tarla), plot_df$treatment_status,
                                   ifelse(plot_df$tarla == "pre", plot_df$treatment_status, "Tarla"))

# plot_df$treatment_status <- ifelse(plot_df$treatment_status == "naive", "Naive", "SOC")

plot_df$treatment_status <- factor(plot_df$treatment_status, levels = c("Naive", "CTX + ICI","Tarla"))

# Remove post-tarla samples
# plot_df <- plot_df %>%
#   filter(tarla != "post" | is.na(tarla))

################################################################################
# Plot violin plots
################################################################################

plot_df$collection_id <- as.character(plot_df$collection_id)


plot_df$color <- ifelse(plot_df$condition == "normal", "Normal", plot_df$collection_id)

p2 <- ggviolin(plot_df, x="collection_id" ,y="expression", fill="color", lwd=.3, outlier.size = .1,draw_quantiles =0.5,)+
  # stat_pvalue_manual(stat.test, y.position = c(6.2,7,6.5),label = "p.adj.signif",size=5,tip.length = 0)+
  stat_compare_means(label = "p.format", size=5)+
  facet_wrap(~antigen,nrow=2)+
  ylim(0,NA)+
  labs(y="Expression",
       x= "")+
  scale_fill_manual(name = "Treatment",values=c("#1f77b3","#fe7e08","#2ba02d","gray"))+
  theme(axis.title = element_text(size=20),
        axis.text.x = element_text(size=16,angle=45,hjust=1,vjust=1),
        axis.text.y = element_text(size=16),
        strip.text = element_text(face = "bold", size=16), 
        strip.background = element_blank(),
        legend.position = "none")

p2


plot_df$collection_id <- ifelse(plot_df$condition == "normal", "Normal", plot_df$collection_id)

stat.test <- plot_df %>%
  wilcox_test(expression ~ collection_id, 
              comparisons = list(c("Normal","MDA-SC547-1"),c("Normal","MDA-SC547-2"),c("Normal","MDA-SC547-3")),
              p.adjust.method = "BH") %>%
  add_significance() 


p3 <- ggviolin(plot_df, x="collection_id" ,y="expression", fill="color", lwd=.3, outlier.size = .1,draw_quantiles =0.5,
              )+
  stat_pvalue_manual(stat.test, y.position = c(6.2,6.5,7),label = "p.adj.signif",size=5,tip.length = 0)+
  facet_wrap(~antigen,nrow=2)+
  ylim(0,NA)+
  labs(y="Expression",
       x= "")+
  scale_color_manual(name = "Treatment",values=c("#1f77b3","#fe7e08","#2ba02d","gray"))+
  scale_fill_manual(name = "Treatment",values=c("#1f77b3","#fe7e08","#2ba02d","gray"))+
  theme_classic()+
  theme(axis.title = element_text(size=20),
        axis.text.x = element_text(size=16,angle=45,hjust=1,vjust=1),
        axis.text.y = element_text(size=16),
        strip.text = element_text(face = "bold", size=16), 
        strip.background = element_blank(),
        legend.position = "none")
p3

p4 <- ggviolin(plot_df, x = "collection_id", y = "expression",
         fill = "color", color = NA, alpha = 0.5,
         trim = F, width = 1) +
  geom_violin(data=plot_df,
              aes(x = collection_id, y = expression), color="black",
              fill = NA, linewidth = .5,
              trim = F, width = 1, position = position_dodge(0.9), draw_quantiles = 0.5)+
  stat_pvalue_manual(stat.test, y.position = c(6.2,6.5,7),label = "p.adj.signif",size=5,tip.length = 0)+
  facet_wrap(~antigen,nrow=2)+
  ylim(0,NA)+
  labs(y="Expression",
       x= "")+
  scale_color_manual(name = "Treatment",values=c("#1f77b3","#fe7e08","#2ba02d","gray"))+
  scale_fill_manual(name = "Treatment",values=c("#1f77b3","#fe7e08","#2ba02d","gray"))+
  theme_classic()+
  theme(axis.title = element_text(size=20),
        axis.text.x = element_text(size=16,angle=45,hjust=1,vjust=1),
        axis.text.y = element_text(size=16),
        strip.text = element_text(face = "bold", size=16), 
        strip.background = element_blank(),
        legend.position = "none")
p4

p5 <- ggviolin(plot_df, x = "collection_id", y = "expression",
               fill = "color", color = NA, alpha = 0.5,
               trim = F, width = 1) +
  geom_violin(data=plot_df,
              aes(x = collection_id, y = expression, color=color),
              fill = NA, linewidth = .5,
              trim = F, width = 1, position = position_dodge(0.9), draw_quantiles = 0.5)+
  stat_pvalue_manual(stat.test, y.position = c(6.2,6.5,7),label = "p.adj.signif",size=5,tip.length = 0)+
  facet_wrap(~antigen,nrow=2)+
  ylim(0,NA)+
  labs(y="Expression",
       x= "")+
  scale_color_manual(name = "Treatment",values=c("#1f77b3","#fe7e08","#2ba02d","gray"))+
  scale_fill_manual(name = "Treatment",values=c("#1f77b3","#fe7e08","#2ba02d","gray"))+
  theme_classic()+
  theme(axis.title = element_text(size=20),
        axis.text.x = element_text(size=16,angle=45,hjust=1,vjust=1),
        axis.text.y = element_text(size=16),
        strip.text = element_text(face = "bold", size=16), 
        strip.background = element_blank(),
        legend.position = "none")

p3

tiff(glue("figures/normal_vs_SC547_DLL3_expression_violinplot_v1.tiff"), width=150,height=165, units = "mm", res=600)
print(p3)
dev.off()


tiff(glue("figures/normal_vs_SC547_DLL3_expression_violinplot_v2.tiff"), width=150,height=165, units = "mm", res=600)
print(p4)
dev.off()


tiff(glue("figures/normal_vs_SC547_DLL3_expression_violinplot_v3.tiff"), width=150,height=165, units = "mm", res=600)
print(p5)
dev.off()
p4
p5
