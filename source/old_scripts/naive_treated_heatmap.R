# Plots heatmap of median protein expression in cluster 1 with treatment status anno
source("source/cytof_de_function.R")


sce <- readRDS("data/cytof_objects/sclc_cytof_sce_object_with_clusters.rds")

sce <- sce[,colData(sce)$condition == "cancer"]

ctcs <- sce[,colData(sce)$new_clusters == 1]

ctcs <- ctcs[rowData(ctcs)$marker_class == "state",]

ctcs <- CATALYST::cluster(ctcs, features = "state",
                          xdim = 10, ydim = 10, maxK = 20, seed = script_seed)


ctcs <- runDR(ctcs, "UMAP", cells = 5e3, features = "state")



ctcs@metadata$delta_area

plotDR(ctcs, "UMAP", color_by = "meta10", scale = T)+
  geom_point(size=1)

plotDR(ctcs, "UMAP", color_by = "treatment_status", scale = T)+
  geom_point(size=1)


plotDR(ctcs, "UMAP", color_by = "POU2F3", scale=T)


colData(ctcs)$new_clusters <- cluster_ids(ctcs, "meta10")
################################################################################



################################################################################
# Differential expression

df <- cytof_de(ctcs, method = "wilcox", metric = "median")

plot_df <- df %>% 
  mutate(significance = factor(ifelse(p_adj < .05, "signif", "ns"), levels=c("signif","ns"))) %>% 
  dplyr::filter(cluster %in% signif_cluster) %>% 
  mutate(protein = reorder_within(protein, logfc, cluster))


plot_df <- df %>% 
  mutate(significance = factor(ifelse(p_adj < .05, "signif", "ns"), levels=c("signif","ns"))) %>% 
  mutate(protein = reorder_within(protein, logfc, cluster))

ggplot(plot_df)+
  geom_col(aes(x=as.numeric(logfc), y=protein, fill=significance))+
  facet_wrap(~cluster, scales = "free_y")+
  scale_y_reordered()+
  ylab("Protein")+
  xlab("Median log2FC")



################################################################################
y <- assay(ctcs, "exprs")

df <- data.frame(t(y), colData(ctcs), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")
gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(ctcs)))

################################################################################
# Patient Heatmap
################################################################################
temp <- gg_df %>% 
  group_by(sample_id, antigen) %>% 
  summarise(median(expression)) %>% 
  pivot_wider(names_from = "antigen", values_from = "median(expression)") %>% 
  column_to_rownames("sample_id")


temp_scaled <- apply(temp, MARGIN = 2, FUN = function(x) scale_values(x))


Heatmap(temp_scaled, name = "Median\nExpression", col = col_fun)

dim(ctcs)

heatmap_metadata <- data.frame(colData(ctcs)) %>% 
  dplyr::select(sample_id,sample_num,treatment_status) %>% 
  distinct()

# Create Annotations
treatment_anno <-  rowAnnotation("Status" = as.vector(heatmap_metadata$treatment_status), 
                                 col = list("Status"=c("treated" = "darkorange","naive" = "royalblue")),
                                 show_annotation_name = T)


Heatmap(temp_scaled, name = "Median\nExpression", col = col_fun, right_annotation = treatment_anno)


temp_scaled_tf <- temp_scaled[,colnames(temp_scaled) %in% c("POU2F3","NeuroD1","ASCL1")]
Heatmap(temp_scaled_tf, name = "Median\nExpression", col = col_fun)



