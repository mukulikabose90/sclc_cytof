source("source/sclc_cytof_functions.R")

set.seed(42)
################################################################################
# Read in data
################################################################################

sce <- readRDS("data/cytof_objects/sclc_all_samples_object_no_qc.rds")

################################################################################
# Remove cell line and PE samples
################################################################################
blood_samples <- as.data.frame(sce@colData) %>%
  dplyr::filter(sample_type == "blood") %>%
  pull(collection_id) %>%
  as.character()

sce <- sce[,sce$collection_id %in% blood_samples]

################################################################################
# Proportional down-sampling
################################################################################
sce$cell_id <- 1:ncol(sce)

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
# Get protein markers
################################################################################
marker_info <- read.csv("data/cytof_panel_info.csv")
marker_info <- data.frame(marker_info, stringsAsFactors = FALSE)

markers <- marker_info %>%
  dplyr::filter(marker_class != "none") %>%
  pull(antigen)

saveRDS(markers, "data/protein_markers.rds")

################################################################################
# cyCombine batch correction
################################################################################
# Set replicates as anchor
replicates <- as.data.frame(colData(sce)) %>% 
  dplyr::count(collection_id,experiment_id) %>% 
  dplyr::count(collection_id) %>% 
  filter(n > 1) %>% 
  pull(collection_id) %>% 
  as.character()

sce$anchor <- ifelse(as.character(sce$collection_id) %in% replicates, as.character(sce$collection_id), as.character(sce$sample_id))

# get expression data
y <- assay(sce, "exprs")

# Set cell IDs
colnames(y) <- paste0("cell_",1:ncol(y))

# Merge metadata and expression
df <- data.frame(t(y), colData(sce), check.names = FALSE)

# Change experiment id column to batch
colnames(df)[which(colnames(df) == "experiment_id")] <- "batch"

# Run batch correction function
corrected <- batch_correct(df,
                           markers = markers,
                           anchor = "anchor")

# Remove metadata columns. Only proteins remain
corrected <- as.matrix(t(corrected[,3:40]))

# Reorder cells to match orignal order
corrected <- corrected[rownames(y),]

# Replace old expression matrix with corrected
assay(sce, "exprs") <- corrected

################################################################################

saveRDS(sce, "data/cytof_objects/sclc_all_samples_object.rds")

