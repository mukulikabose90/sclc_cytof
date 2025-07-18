source("source/sclc_cytof_functions.R")

set.seed(42)
################################################################################
# Read in data
################################################################################

sce <- readRDS("data/cytof_objects/sclc_all_samples_object_no_qc.rds")

################################################################################
# Get protein markers
################################################################################
marker_info <- read.csv("data/cytof_panel_info.csv")
marker_info <- data.frame(marker_info, stringsAsFactors = FALSE)

state_markers <- marker_info %>%
  dplyr::filter(marker_class == "state") %>%
  pull(antigen)

saveRDS(state_markers, "data/state_markers.rds")

################################################################################
# Remove outlier experiment
################################################################################
sce <- sce[,sce$experiment_id != "531050"]

################################################################################
# Remove low viability samples
################################################################################
low_viability_samples <- paste0("NORMAL", c(10,11,15,18))

sce <- sce[,!sce$patient_id %in% low_viability_samples]

################################################################################
# Remove non-SCLC samples
################################################################################
NSCLC_samples <- c("SC378","SC368","SC404","SC383")

sce <- sce[,!sce$patient_id %in% NSCLC_samples]
################################################################################
# cyCombine batch correction
################################################################################
# markers_to_use <- state_markers
# 
# # Set replicates as anchor
# replicates <- as.data.frame(colData(sce)) %>%
#   dplyr::count(collection_id,experiment_id) %>%
#   dplyr::count(collection_id) %>%
#   dplyr::filter(n > 1) %>%
#   pull(collection_id) %>%
#   as.character()
# 
# sce$anchor <- ifelse(as.character(sce$collection_id) %in% replicates, as.character(sce$collection_id), as.character(sce$sample_id))
# 
# # get expression data
# y <- assay(sce, "exprs")
# 
# # Set cell IDs
# colnames(y) <- paste0("cell_",1:ncol(y))
# 
# y <- y[rownames(y) %in% markers_to_use,]
# 
# # Merge metadata and expression
# df1 <- data.frame(t(y), colData(sce)[which(colnames(colData(sce)) %in% c("collection_id","experiment_id"))], check.names = FALSE)
# 
# #df1
# colnames(df1)[which(colnames(df1) == "experiment_id")] <- "batch"
# colnames(df1)[which(colnames(df1) == "collection_id")] <- "sample"
# 
# head(t(y))
# 
# df1$sample <- as.character(df1$sample)
# df1$sample <- sapply(df1$sample, FUN = function(x) strsplit(x, "-")[[1]][1])
# df1$sample <- gsub("SC","",df1$sample)
# df1$sample <- gsub("NORMAL","",df1$sample)
# 
# df1$sample <- gsub("[A-Za-z]", "", df1$sample)
# 
# 
# 
# df1$batch <- as.numeric(df1$batch)
# df1$sample <- as.numeric(df1$sample)
# 
# colnames(df1) <- gsub("-","",colnames(df1))
# colnames(df1) <- gsub("_","",colnames(df1))
# 
# markers <- gsub("-","",markers_to_use)
# markers <- gsub("_","",markers)
# 
# # detect_batch_effect(df1,
# #                     batch_col = 'batch',
# #                     out_dir = paste0("figures/batch_effects/"),
# #                     seed = 434,
# #                     name = 'CyTOF Data',
# #                     markers = markers)
# # 
# # detect_batch_effect_express(df1, downsample = 10000, out_dir = 'figures/batch_effects')
# 
# 
# 
# ################################################################################
# df <- data.frame(t(y), colData(sce)[which(colnames(colData(sce)) %in% c("collection_id","experiment_id","condition"))], check.names = FALSE)
# 
# 
# 
# # Change experiment id column to batch
# colnames(df)[which(colnames(df) == "experiment_id")] <- "batch"
# colnames(df)[which(colnames(df) == "collection_id")] <- "sample"
# 
# head(t(y))
# 
# df$sample <- as.character(df$sample)
# df$sample <- sapply(df$sample, FUN = function(x) strsplit(x, "-")[[1]][1])
# df$sample <- gsub("SC","",df$sample)
# df$sample <- gsub("NORMAL","",df$sample)
# 
# 
# df$batch <- as.numeric(df$batch)
# df$sample <- as.numeric(df$sample)
# 
# colnames(df) <- gsub("-","",colnames(df))
# colnames(df) <- gsub("_","",colnames(df))
# 
# 
# # markers <- c("NeuroD1", "PDL1", "DLL3", "pYAP", "CD44", "ECad", "CD24", "POU2F3", "EpCAM", "ASCL1")
# 
# # Run batch correction function
# corrected <- batch_correct(df,
#                            markers = markers,
#                            covar = "condition")
# 
# # Remove metadata columns. Only proteins remain
# df1 <- corrected[,3:ncol(corrected)]
# corrected <- corrected[,3:18]
# colnames(corrected) <- markers_to_use
# 
# 
# dim(corrected)
# # corrected <- corrected[,3:ncol(corrected)]
# 
# 
# # detect_batch_effect(df1,
# #                     batch_col = 'batch',
# #                     out_dir = paste0("figures/post_batch_effects/"),
# #                     seed = 434,
# #                     name = 'CyTOF Data',
# #                     markers = markers)
# # 
# # detect_batch_effect_express(df1, downsample = 10000, out_dir = 'figures/post_batch_effects')
# 
# 
# assay(sce, "exprs")[colnames(corrected),] <- t(corrected)

################################################################################
# Proportional down-sampling
################################################################################
sce$cell_id <- paste0("cell_",1:ncol(sce))

original_sce <- sce

all_cell_ids <- sce$cell_id

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

remaining_cell_ids <- sce$cell_id

left_out_cell_ids <- all_cell_ids[!all_cell_ids %in% remaining_cell_ids]

removed_sce <- original_sce[,original_sce$cell_id %in% left_out_cell_ids]

saveRDS(removed_sce, "data/cytof_objects/sce_not_sampled.rds")

################################################################################

saveRDS(sce, "data/cytof_objects/sclc_all_samples_object.rds")

