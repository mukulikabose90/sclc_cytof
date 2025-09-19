source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################

sce <- readRDS("data/cytof_objects/all_samples_ctcs_with_subtype.rds")
sce <- readRDS("data/cytof_objects/tarla_sce.rds")

colData(sce)$condition <- factor(colData(sce)$condition, levels=c("normal", "cancer"))
sce@metadata$experiment_info$condition <- factor(sce@metadata$experiment_info$condition, levels=c("normal", "cancer"))

sce <- sce[,colData(sce)$condition == "cancer"]
sce <- sce[,!is.na(colData(sce)$tarla)]

dim(sce)

as.data.frame(colData(sce)) %>% 
  dplyr::count(tarla,collection_id) %>% 
  distinct() %>% 
  dplyr::count(tarla)



sce <- CATALYST::cluster(sce, features = "state",
                         xdim = 10, ydim = 10, maxK = 20, seed = script_seed)


sce <- runDR(sce, "UMAP", cells = 5e3, features = "state")

markers_to_use <- rowData(sce)$marker_name[rowData(sce)$marker_class == "state"]
y <- assay(sce, "exprs")
y <- t(y[markers_to_use,])
y <- as.data.frame(apply(y, 2, function(x) (x - min(x)) / (max(x) - min(x))))

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
  
  
  curr_plot
  
  jpeg(glue("figures/marker_expression_umaps/tarla_{curr_marker}.jpg"), width=120,height=100, units = "mm", res=1000)
  print(curr_plot)
  dev.off()
}
