source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################
# Read in CyTOF data with cluster assignments
################################################################################
sce <- readRDS("data/cytof_objects/sclc_all_samples_with_clusters.rds")

cancer_enriched_clusters <- readRDS("data/cancer_enriched_clusters.rds")

# Subset to cancer cells in cancer_enriched cluster
cancer_enriched <- sce[,colData(sce)$new_clusters %in% cancer_enriched_clusters]
cancer_enriched <- cancer_enriched[,colData(cancer_enriched)$condition == "cancer"]

# Subset to only state markers
cancer_enriched <- cancer_enriched[rowData(cancer_enriched)$marker_class == "state",]

# Give each cell an ID
colData(cancer_enriched)$cell_id <- 1:nrow(colData(cancer_enriched))
colData(cancer_enriched)$cell_id <- paste0("cell_",1:nrow(colData(cancer_enriched)))

################################################################################
# Cluster cells using FlowSOM then run UMAP
################################################################################
markers_to_use <- readRDS("data/state_markers.rds")
# markers_to_use <- markers_to_use[-which(markers_to_use %in% c("NeuroD1","ASCL1","POU2F3","p-YAP","SLUG","Twist","p-Rb"))]

cancer_enriched <- CATALYST::cluster(cancer_enriched, features = markers_to_use,
                                     xdim = 10, ydim = 10, maxK = 20, seed = script_seed)

cancer_enriched <- runDR(cancer_enriched, "UMAP", cells = 5e3, features = markers_to_use)


markers_to_use <- rowData(cancer_enriched)$marker_name[rowData(cancer_enriched)$marker_class == "state"]
y <- assay(cancer_enriched, "exprs")
y <- t(y[markers_to_use,])
# y <- as.data.frame(apply(y, 2, function(x) (x - min(x)) / (max(x) - min(x))))

xy <- reducedDim(cancer_enriched, "UMAP")[, 1:2]
colnames(xy) <- c("x", "y")

df <- data.frame(colData(cancer_enriched), xy,y, check.names = FALSE)


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
  
  jpeg(glue("figures/marker_expression_umaps/ce_clusters_{curr_marker}_expression.jpg"), width=120,height=100, units = "mm", res=1000)
  print(curr_plot)
  dev.off()
}

