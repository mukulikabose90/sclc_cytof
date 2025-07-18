source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################
sce <- readRDS("data/cytof_objects/sclc_all_samples_with_clusters.rds")

cancer_enriched_clusters <- readRDS("data/cancer_enriched_clusters.rds")

# sce$cancer_enriched <- ifelse(sce$new_clusters %in% cancer_enriched_clusters & sce$condition == "cancer", "cancer_enriched", "non-cancer_enriched")
sce$cancer_enriched <- ifelse(sce$new_clusters %in% cancer_enriched_clusters, "cancer_enriched", "non-cancer_enriched")

markers_to_use <- c("NeuroD1","ASCL1","POU2F3","p-Rb")



p1 <- create_marker_boxplots(sce, markers_to_use, "new_clusters","cancer_enriched")

p1 <- p1+
  labs(y="Expression",
       x="",
       fill="")+
  scale_fill_manual(
    values = c("cancer_enriched" = "#E57373", "non-cancer_enriched" = "#64B5F6"),
    labels = c("cancer_enriched" = "Cancer Enriched", "non-cancer_enriched" = "Non-Cancer Enriched"))+
  rremove("legend")




p2 <- create_marker_boxplots(sce, markers_to_use, "cancer_enriched","cancer_enriched")

p2 <- p2+
  stat_compare_means(label = "p.signif", tip.length = 0, comparisons = list(c("cancer_enriched","non-cancer_enriched")))+
  labs(y="Expression",
       x="",
       fill="")+
  ylim(0,10)+
  scale_fill_manual(
    values = c("cancer_enriched" = "#E57373", "non-cancer_enriched" = "#64B5F6"),
    labels = c("cancer_enriched" = "Cancer Enriched", "non-cancer_enriched" = "Non-Cancer Enriched"))+
  scale_x_discrete(labels = c("cancer_enriched" = "Cancer Enriched", "non-cancer_enriched" = "Normal Enriched"))+
  rremove("legend")

p2

markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam", "E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist", "SLUG", "PD-L1", "p-YAP", "CD44", "CD24")

p3 <- create_marker_boxplots(sce, markers_to_use, "new_clusters","cancer_enriched")

p3 <- p3+
  labs(y="Expression",
       x="",
       fill="")+
  scale_fill_manual(
    values = c("cancer_enriched" = "#E57373", "non-cancer_enriched" = "#64B5F6"),
    labels = c("cancer_enriched" = "Cancer Enriched", "non-cancer_enriched" = "Non-Cancer Enriched"))+
  rremove("legend")

################################################################################

jpeg(glue("figures/all_samples_cluster_tf_expression_boxplot.jpg"), width=160,height=160, units = "mm", res=1000)
print(p1)
dev.off()

jpeg(glue("figures/all_samples_ce_vs_nonce_tf_expression_boxplot.jpg"), width=160,height=160, units = "mm", res=1000)
print(p2)
dev.off()

jpeg(glue("figures/all_samples_cluster_expression_boxplot.jpg"), width=320,height=200, units = "mm", res=1000)
print(p3)
dev.off()
