source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################
# Read in CyTOF data with cluster assignments
################################################################################
sce <- readRDS("data/cytof_objects/sclc_all_samples_with_clusters.rds")

# sce <- sce[,sce$patient_id != "SC443"]


ctc_clusters <- readRDS("data/ctc_clusters.rds")
# ctc_clusters <- c(6)

# Subset to cancer cells in CTC cluster
ctcs <- sce[,colData(sce)$new_clusters %in% ctc_clusters]
ctcs <- ctcs[,colData(ctcs)$condition == "cancer"]

# Subset to only state markers
ctcs <- ctcs[rowData(ctcs)$marker_class == "state",]

# Give each cell an ID
colData(ctcs)$cell_id <- 1:nrow(colData(ctcs))
colData(ctcs)$cell_id <- paste0("cell_",1:nrow(colData(ctcs)))




################################################################################
# Cluster cells using FlowSOM then run UMAP
################################################################################
markers_to_use <- readRDS("data/state_markers.rds")
# markers_to_use <- markers_to_use[-which(markers_to_use %in% c("NeuroD1","ASCL1","POU2F3","p-YAP","SLUG","Twist","p-Rb"))]

ctcs <- CATALYST::cluster(ctcs, features = markers_to_use,
                         xdim = 10, ydim = 10, maxK = 20, seed = script_seed)

ctcs <- runDR(ctcs, "UMAP", cells = 5e3, features = markers_to_use)

################################################################################
# Identify optimal number of clusters and assign cells
################################################################################
ctcs@metadata$delta_area


CATALYST::plotDR(ctcs, color_by = "meta6")

CATALYST::plotDR(ctcs, color_by = "meta6",facet_by = "patient_id")

CATALYST::plotDR(ctcs, color_by = "patient_id")


colData(ctcs) %>% 
  as.data.frame() %>% 
  dplyr::count(patient_id, new_clusters)


colData(ctcs)$new_clusters <- cluster_ids(ctcs, "meta9")


y <- assay(ctcs, "exprs")

# 
# fviz_nbclust(y, kmeans, method='silhouette')+
#   ggtitle("Optimal Number of Subclusters (CTCs)")
# 
# 
# 
# metric_to_use <- "median"
# 
# df <- cytof_de(ctcs, method = "wilcox", metric = metric_to_use, ident = "new_clusters")
# 
# ################################################################################
# # Create DE barplots
# ################################################################################
# plot_df <- df %>% 
#   mutate(significance = factor(ifelse(p_adj < .05, "*", ""), levels=c("*",""))) %>% 
#   mutate(protein = reorder_within(protein, logfc, ident_list)) %>% 
#   mutate(star_x = ifelse(logfc > 0, logfc+.25,logfc-.5)) %>% 
#   mutate(up_down= ifelse(logfc>0,"up","down"))
# 
# # Reorder cluster factor
# plot_df$ident_list <- paste0("Cluster ",plot_df$ident_list)
# plot_df$ident_list <- factor(plot_df$ident_list, levels=paste0("Cluster ", c(1:10)))
# 
# # Dataframe for plotting only CTC clusters
# ctc_df <- plot_df %>% 
#   dplyr::filter(ident_list %in% paste0("Cluster ", ctc_clusters))
# 
# 
# x_axis_label <- gsub("m","M",metric_to_use)
# 
# # Create plot for CTC clusters
# p1 <- ggplot(plot_df,aes(x=as.numeric(logfc), y=protein, fill=logfc))+
#   geom_col(color="darkgray",linewidth=.001)+
#   geom_text(aes(x=star_x, label=significance), size = 3)+
#   facet_wrap(~ident_list, scales = "free_y",nrow=1)+
#   scale_y_reordered()+
#   labs(fill = "")+
#   guides(fill="none")+
#   # scale_fill_manual(values=c("royalblue4", "firebrick"))+
#   scale_fill_gradient2(low = "royalblue4", mid = "white", high = "firebrick", midpoint = 0)+
#   ylab("Protein")+
#   xlab(glue("{x_axis_label} log(FC)"))+
#   theme_classic()+
#   theme(panel.grid.minor = element_blank(), 
#         strip.text = element_text(face = "bold", size=8), 
#         axis.text = element_text(color = "black", size=6),
#         axis.title = element_text(size=8),
#         legend.text = element_text(size=6),
#         legend.title = element_text(size=8))
# 
# p1



df <- data.frame(t(y), colData(ctcs), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")

gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(sce)))

sclc_tfs <- c("NeuroD1","ASCL1","POU2F3","p-Rb")

sclc_tfs <- markers

temp <- gg_df %>%
  dplyr::filter(antigen %in% sclc_tfs)




p <- ggboxplot(temp, x="new_clusters",y="expression", fill="new_clusters")

  + facet_wrap(~antigen, scales = "free_y")


sum(ctcs$new_clusters == 1)
  

ctcs <- ctcs[,ctcs$new_clusters %in% c(2,3,4,5)]



# Select expression data
y <- assay(ctcs, "exprs")

#Create tidy dataframe for each sample
df <- data.frame(t(y), colData(ctcs), check.names = FALSE)
value <- ifelse("exprs" == "exprs", "expression", "exprs")
gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(ctcs)))

# Create heatmap of protein expression for each cell
heatmap <- gg_df %>% 
  select(cell_id,antigen,expression) %>% 
  pivot_wider(names_from = antigen, values_from = expression) %>% 
  column_to_rownames("cell_id")

# Standardize expression within each protein
heatmap_scaled <- scale(heatmap)

# Subset to only SCLC subtype TFs 
heatmap_tf_scaled <- heatmap_scaled[,colnames(heatmap_scaled) %in% c("POU2F3","NeuroD1","ASCL1")]

#Set up metadata table to create annotations
heatmap_metadata <- data.frame(colData(ctcs)) %>% 
  dplyr::select(cell_id,collection_id,sample_id,patient_id,condition,sample_num, treatment_status) %>% 
  select(-c(sample_id)) %>% 
  distinct()

all_samples_heatmap <- t(heatmap_tf_scaled)

#############################################################################
# Create heatmap with all samples
#############################################################################
col_fun = colorRamp2(c(-2, 0, 2), c("royalblue4","lightgray", "firebrick4"))

# Create collection ID annotation
colors_to_use <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(heatmap_metadata$collection_id)))
names(colors_to_use) <- unique(heatmap_metadata$collection_id)
sample_anno <- HeatmapAnnotation("Sample ID" = heatmap_metadata$collection_id, 
                                 col = list("Sample ID"= colors_to_use),
                                 show_annotation_name = T,
                                 annotation_legend_param = list(ncol=3))

# Create patient ID annotation
colors_to_use <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(heatmap_metadata$patient_id)))
names(colors_to_use) <- unique(heatmap_metadata$patient_id)
patient_anno <- HeatmapAnnotation("Patient ID" = heatmap_metadata$patient_id, 
                                  col = list("Patient ID"= colors_to_use),
                                  show_annotation_name = T,
                                  annotation_legend_param = list(ncol=2))


ht <- Heatmap(all_samples_heatmap, column_km = 4, top_annotation = sample_anno, name="Expression",
              cluster_columns = T, cluster_rows = F, show_column_names=F,col = col_fun,
              row_dend_reorder = F, column_title = c("SCLC-A", "SCLC-P","SCLC-I","SCLC-N"))

all_samples_ht <- draw(ht)





fviz_nbclust(t(all_samples_heatmap), kmeans, method='silhouette')+
  ggtitle("Optimal Number of Subclusters (CTCs)")


clusters <- column_order(all_samples_ht)

# Rename clusters
names(clusters) <- c("A","P","I","N")
# names(clusters) <- c("P","A","I","N")

# Create dataframe of cell IDs and associated subtypes
subtypes_df <- list()
for(curr_subtype in names(clusters)){
  subtypes_df <- append(subtypes_df, list(cbind(colnames(all_samples_heatmap[,clusters[[curr_subtype]]]),curr_subtype)))
}
subtypes_df <- do.call(rbind,subtypes_df)
colnames(subtypes_df) <- c("cell_id","subtype")

# Add original order of cell in colData, then merge based on cell ID. Finally rearragne to orignal order
colData(ctcs)$original_order  <- 1:nrow(colData(ctcs))
subtype_order <- merge(colData(ctcs), subtypes_df, by="cell_id")
subtype_order <- subtype_order[order(subtype_order$original_order), ]

# Check that cell order is back to orginal order
all(subtype_order$cell_id == colData(ctcs)$cell_id)

# Add subtypes to colData
colData(ctcs)$subtype <- subtype_order$subtype

# Save data with subtype assignments
saveRDS(ctcs, "data/cytof_objects/all_samples_ctcs_with_subtype.rds")

##################################

markers_to_use <- rowData(sce)$marker_name[rowData(sce)$marker_class == "state"]
y <- assay(sce, "exprs")
y <- t(y[markers_to_use,])
# y <- as.data.frame(apply(y, 2, function(x) (x - min(x)) / (max(x) - min(x))))

xy <- reducedDim(sce, "UMAP")[, 1:2]
colnames(xy) <- c("x", "y")

df <- data.frame(colData(sce), xy,y, check.names = FALSE)


markers <- colnames(y)

for(curr_marker in markers){
  cat(curr_marker,"\n")
  
  # Plot UMAP
  curr_plot <- ggplot(df)+
    geom_point(aes(x=x, y=y, color=!!sym(curr_marker)),size=.01)+
    xlab("UMAP 1")+
    ylab("UMAP 2")+
    labs(color = curr_marker)+
    scale_color_gradientn(colors = c("lightblue","red2"))+
    theme_classic() +
    theme(panel.grid.minor = element_blank(), 
          strip.text = element_text(face = "bold", size=8), 
          axis.text = element_text(color = "black", size=8),
          axis.title = element_text(size=8),
          legend.text = element_text(size=6),
          legend.title = element_text(size=8))
  
  
  # curr_plot
  
  jpeg(glue("figures/marker_expression_umaps/ctcs_{curr_marker}.jpg"), width=120,height=100, units = "mm", res=1000)
  print(curr_plot)
  dev.off()
}


