

ctcs <- readRDS("data/cytof_objects/all_samples_ctcs_with_subtype.rds")

subtypes <- c("A","N","P","I")
cols <- list()
for(curr_subtype in subtypes){
  curr_sce <- ctcs[,ctcs$subtype == curr_subtype]
  
  
  pre <- rowMeans(curr_sce@assays@data$exprs[,curr_sce$treatment_status == "naive"])
  post <- rowMeans(curr_sce@assays@data$exprs[,curr_sce$treatment_status == "treated"])
  
  # cols <- append(cols,list((post+.000001)/(pre+.000001)))
  cols <- append(cols,list((post)/(pre)))
  
}

heatmap <- do.call(cbind,cols)

colnames(heatmap) <- subtypes
rownames(heatmap) <- rownames(ctcs)


Heatmap(log2(heatmap),cluster_rows = F,cluster_columns = F, name="log2(FC)",
        column_names_rot = 0)


######################################################


################################################################################
# Read in data
################################################################################
ctcs <- readRDS("data/cytof_objects/all_samples_ctcs_with_subtype.rds")

subtypes <- c("A","N","P","I")
heatmaps <- list()
for(curr_subtype in subtypes){
  curr_sce <- ctcs[,ctcs$subtype == curr_subtype]
  
  pre <- rowMeans(curr_sce@assays@data$exprs[,curr_sce$treatment_status == "naive"])
  post <- rowMeans(curr_sce@assays@data$exprs[,curr_sce$treatment_status == "treated"])
  
  # cols <- append(cols,list((post+.000001)/(pre+.000001)))
  curr_heatmap <- cbind((pre)/(pre),(post)/(pre))
  
  col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  ht <- Heatmap(log2(curr_heatmap),cluster_rows = F,cluster_columns = F, column_title = curr_subtype,
                col = col_fun, name = "log2(FC)",rect_gp = gpar(col = "gray", lwd = 2))
  
  heatmaps <- append(heatmaps, list(ht))
}

heatmaps[[1]]+heatmaps[[2]]+heatmaps[[3]]+heatmaps[[4]]









