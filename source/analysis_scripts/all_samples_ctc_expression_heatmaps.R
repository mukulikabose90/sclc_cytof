source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################

ctcs <- readRDS("data/cytof_objects/all_samples_ctcs_with_subtype.rds")

markers_to_use <- rowData(ctcs)$marker_name[rowData(ctcs)$marker_class == "state"]

# markers_to_use <- c("E-Cad","Vimentin","CD44","CD24","Twist","MUC-1")


# Subset markers
ctcs <- ctcs[markers_to_use,]

clusters <- as.numeric(unique(ctcs$new_clusters))

y <- assay(ctcs, "exprs")
dim(t(y))


y <- scale(t(y))

dim(y)



ht_list <- list()
for(curr_cluster in clusters){
  
  curr_mat <- y[ctcs$new_clusters == curr_cluster,]
  
  # scaled_matrix <- scale(y)
  
  
  
  col_fun = colorRamp2(c(min(curr_mat), 0.5, max(curr_mat)), c("royalblue", "white", "firebrick2"))
  
  ht <- Heatmap(curr_mat, cluster_rows = F,cluster_columns = F,
                column_title = paste0("Cluster ", curr_cluster),
                name=" ",show_heatmap_legend = FALSE, col = col_fun)
  ht
  
  ht_list <- append(ht_list, list(ht))
  
}


at = seq(0, 1, by = 0.2)

col_fun = colorRamp2(c(0, 0.5, 1), c("firebrick", "white", "royalblue4"))
lgd = Legend(at = rep("",length(at)), title = "", legend_gp = gpar(fill = col_fun(at)))



draw(lgd)



grid.newpage()

# Create viewports
pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 3)))

# Draw heatmap 1
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(ht_list[[1]], newpage = FALSE)
upViewport()

# Draw heatmap 2
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(ht_list[[2]], newpage = FALSE)
upViewport()



# Draw legend
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
draw(lgd)
upViewport()
