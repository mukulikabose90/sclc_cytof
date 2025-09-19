source("source/cytof_de_function.R")
script_seed <- 42
set.seed(script_seed)


sce <- readRDS("data/cytof_objects/sclc_cytof_sce_object_with_clusters.rds")

ctcs <- sce[,colData(sce)$new_clusters == 1]

ctcs <- ctcs[rowData(ctcs)$marker_class == "state",]

all_patients <- levels(colData(sce)$patient_id)

all_patients <- all_patients[grepl("SC",all_patients)]

all_patients <- all_patients[which(all_patients != "SC355" & all_patients != "SC370")]

final_heatmap <- list()
for(curr_patient in all_patients){
  curr_sce <- sce[,colData(ctcs)$patient_id == curr_patient]
  
  curr_sce <- CATALYST::cluster(curr_sce, features = "state",
                           xdim = 10, ydim = 10, maxK = 20, seed = script_seed)


  curr_sce <- runDR(curr_sce, "UMAP", cells = 5e3, features = "state")

  # curr_sce@metadata$delta_area
  
  curr_diff <- 10
  values <- curr_sce@metadata$delta_area$data$y
  
  for(i in 2:length(values)){
    prev_val <- values[i-1]
    if(!curr_diff < .01){
      curr_diff <- abs(prev_val-values[i])
      
      # cat(curr_diff,"\n")
      final_i <- i
    }
  }

  # final_i
  
  # plotDR(curr_sce, "UMAP", color_by = glue("meta{final_i}"), scale = T)+
  #   geom_point(size=1)
  
  
  colData(curr_sce)$new_clusters <- cluster_ids(curr_sce, glue("meta{final_i}"))
  
  y <- assay(curr_sce, "exprs")
  
  #Create tidy dataframe for each sample
  df <- data.frame(t(y), colData(curr_sce), check.names = FALSE)
  value <- ifelse("exprs" == "exprs", "expression", "exprs")
  gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
                id.vars = names(colData(curr_sce)))
  
  # Create heatmap of median expression for each patient
  heatmap <- gg_df %>% 
    group_by(new_clusters,antigen) %>% 
    summarise(median(expression)) %>% 
    pivot_wider(names_from = "antigen", values_from = "median(expression)") %>% 
    column_to_rownames("new_clusters")
  
  
  rownames(heatmap) <- paste0(curr_patient, "_", rownames(heatmap))

  final_heatmap <- append(final_heatmap,list(heatmap))
  
}



final_heatmap <- do.call(final_heatmap, rbind)







