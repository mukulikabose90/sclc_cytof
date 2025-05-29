source("source/sclc_cytof_functions.R")

set.seed(42)
################################################################################

sce <- readRDS("data/cytof_objects/sclc_all_samples_object_no_qc.rds")

################################################################################
# Remove cell line samples
################################################################################
blood_samples <- as.data.frame(sce@colData) %>%
  dplyr::filter(sample_type == "blood") %>%
  pull(collection_id) %>%
  as.character()

sce <- sce[,sce$collection_id %in% blood_samples]

################################################################################
# Proportional downsampling
################################################################################
downsample_df <- as.data.frame(colData(sce)) %>% 
  dplyr::count(collection_id,experiment_id) %>% 
  dplyr::count(collection_id) %>% 
  arrange(desc(n))

downsampled_objs <- list()
for(i in 1:nrow(downsample_df)){
  
  curr_sce <- sce[,sce$collection_id == downsample_df[i,1]]
  
  downsample_frac <- 1/downsample_df[i,2]
  
  max_n <- downsample_frac*ncol(curr_sce)
  
  curr_sce <- downsampleSCE(curr_sce,
                            group_by = "collection_id",
                            maxN = max_n)
  
  downsampled_objs <- append(downsampled_objs,list(curr_sce))
}

sce <- do.call(cbind, downsampled_objs)

################################################################################
# markers <- as.data.frame(rowData(sce)) %>%
#   dplyr::filter(marker_class == "state") %>%
#   pull(marker_name)
# 
# temp <- sce[markers,sce$collection_id == "H1105-1"]
# 
# p1 <- plotExprs(temp, color_by = "experiment_id", assay = "exprs")
# 
# p1

# temp <- sce[markers,sce$collection_id == "NJH29-1"]
# 
# p1 <- plotExprs(temp, color_by = "experiment_id", assay = "exprs")
# # 
# temp <- sce_corrected[markers,sce_corrected$collection_id == "NJH29-1"]
# 
# p2 <- plotExprs(temp, color_by = "experiment_id", assay = "exprs")
# 
# p1+ggtitle("NJH29 (no batch correction)")
# p2+ggtitle("NJH29 (batch corrected)")
# 
# ################################################################################

# Get samples that are run in multiple experiments
# to_test <- as.data.frame(sce@colData) %>%
#   dplyr::select(experiment_id,collection_id) %>%
#   distinct() %>%
# dplyr::count(collection_id) %>%
#   arrange(desc(n)) %>%
#   dplyr::filter(n > 1) %>%
#   pull(collection_id) %>%
#   as.character()
# 
# 
# sort(to_test)
# 
# # Checking SC454-1
# temp <- sce[markers,sce$collection_id == "SC454-1"]
# 
# plotExprs(temp, color_by = "experiment_id", assay = "exprs")
# 
# temp <- sce_corrected[markers,sce_corrected$collection_id == "SC414-1"]
# plotExprs(temp, color_by = "experiment_id",assay = "exprs")

################################################################################
# Remove outlier experiments
################################################################################
sce <- sce[,sce$experiment_id != "531050"]
# sce <- sce[,sce$experiment_id != "508814"]
# sce <- sce[,sce$experiment_id != "513549"]

################################################################################
# Remove low viability samples
################################################################################
low_viability_samples <- paste0("NORMAL", c(10,11,15,18))

sce <- sce[,!sce$patient_id %in% low_viability_samples]

################################################################################
# Remove blood bank samples
################################################################################
# blood_bank_samples <- paste0("NORMAL", 7:20)
# 
# length(unique(sce$patient_id))
# 
# sce <- sce[,!sce$patient_id %in% blood_bank_samples]
# 
# as.character(unique(colData(sce)$patient_id))

################################################################################
# remove collections with < 30 cells
################################################################################
samples_to_remove <- as.data.frame(sce@colData) %>%
  dplyr::count(collection_id) %>%
  dplyr::filter(n<30) %>%
  pull(collection_id) %>%
  as.character()

sce <- sce[,!sce$collection_id %in% samples_to_remove]

################################################################################
# limma batch correction
################################################################################

# batch <- as.factor(colData(sce)$experiment_id)
# 
# design <- model.matrix(~ 0 + factor(colData(sce)$condition))  # one-hot encoding of conditions
# colnames(design) <- levels(factor(colData(sce)$condition))
# 
# corrected_exprs <- removeBatchEffect(assay(sce, "exprs"), batch = batch, design = design)
# 
# assay(sce, "exprs") <- corrected_exprs

################################################################################
# cyCombine batch correction
################################################################################
marker_info <- read.csv("data/cytof_panel_info.csv")
marker_info <- data.frame(marker_info, stringsAsFactors = FALSE)

markers <- marker_info %>%
  dplyr::filter(marker_class != "none") %>%
  pull(antigen)

y <- assay(sce, "exprs")

colnames(y) <- paste0("cell_",1:ncol(y))

df <- data.frame(t(y), colData(sce), check.names = FALSE)

colnames(df)[which(colnames(df) == "experiment_id")] <- "batch"

corrected <- batch_correct(df,
                           markers = markers)

corrected <- as.matrix(t(corrected[,3:40]))

corrected <- corrected[rownames(y),]

assay(sce, "exprs") <- corrected

################################################################################

saveRDS(sce, "data/cytof_objects/sclc_all_samples_object.rds")

