source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)

################################################################################

sce <- readRDS("data/cytof_objects/sclc_all_samples_with_clusters.rds")

all_clusters <- as.numeric(sort(unique(sce$new_clusters)))
all_markers <- readRDS("data/state_markers.rds")

heatmap <- matrix(NA, ncol=length(all_markers),nrow=length(all_clusters))
for(i in 1:length(all_clusters)){
  curr_cluster <- all_clusters[i]
  
  for(j in 1:length(all_markers)){
    curr_marker <- all_markers[j]
    
    heatmap[i,j] <- mean(sce@assays@data$exprs[curr_marker,sce$new_clusters == curr_cluster])
  }
}

colnames(heatmap) <- all_markers
rownames(heatmap) <- all_clusters

ht <- Heatmap(scale(heatmap))

draw(ht)