library(ComplexHeatmap)

sce <- readRDS("data/cytof_objects/sclc_cytof_sce_object.rds")
sce <- sce[,colData(sce)$sample_type != "cell_line"]

colData(sce)$cell_num <- 1:ncol(sce)

data <- sce@assays@data$exprs

colnames(data) <- colData(sce)$cell_num

# data.frame(colData(sce)) %>% 
#   dplyr::select(cell_num,patient_id) %>% 
#   deframe()

ha <- HeatmapAnnotation(patient = colData(sce)$patient_id)

Heatmap(cor(data), top_annotation = ha)
