source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)

create_metadata <- function(files,experiment_id){
  
  # This function creates metadata for each file that acts as linker to full metadata csv
  
  filenames <- as.character(files)
  
  # sample information
  data_id <- gsub("_norm_CD45-CD56\\+", "", gsub("\\.fcs$", "", filenames))
  
  data_id <- sapply(data_id, FUN = function(x) strsplit(x,"/")[[1]][4])
  
  data_id <- sapply(data_id, FUN = function(x) strsplit(x,"_")[[1]][1])
  
  data_id <- paste0(data_id,"_",experiment_id)
  
  experiment_metadata <- data.frame(cbind(filenames,data_id))
  rownames(experiment_metadata) <- NULL
  
  return(experiment_metadata)
}

################################################################################
# Read in experiment IDs
################################################################################

all_experiments <- list.files("data/cytof_raw_data")
all_experiments <- all_experiments[-which(all_experiments == "reference_files")]
all_experiments <- all_experiments[!grepl("old",all_experiments)]

################################################################################
# Meta-data: marker information
################################################################################

marker_info <- read.csv("data/cytof_panel_info.csv")
marker_info <- data.frame(marker_info, stringsAsFactors = FALSE)

markers_to_use <- marker_info$fcs_colname

################################################################################
# Read in all CyTOF files
################################################################################

all_metadata <- list()
all_files <- c()
all_flowsets <- list()
for(curr_experiment in all_experiments){
  curr_experiment_id <- strsplit(curr_experiment, "_")[[1]][2]
  
  # Get fcs file names of current experiment
  curr_files <- list.files(glue("data/cytof_raw_data/{curr_experiment}"), pattern = "\\.fcs$", full.names = T)
  
  curr_flowSet <- read.flowSet(curr_files, transformation = FALSE, truncate_max_range = FALSE)
  
  # Subset flowSet to include only selected markers
  # Apply subsetting to each flowFrame, use same markers for each experiment
  curr_flowSet_subset <- fsApply(curr_flowSet, function(ff) {
    ff[, markers_to_use]
  })
  
  # Append current flowset
  all_flowsets <- append(all_flowsets,curr_flowSet_subset)
  
  # Append current file names
  all_files <- append(all_files, curr_files)
  
  # Create and append current metadata
  all_metadata <- append(all_metadata, list(create_metadata(curr_files,curr_experiment_id)))
}

# Merge all flowsets together
final_flowset <- Reduce(rbind2, all_flowsets)

# add clinical metadata
all_metadata <- do.call(rbind, all_metadata)
clinical_metadata <- read.csv("data/cytof_metadata.csv")

# Select the clinical metadata rows and cols you want to keep
clinical_metadata <- clinical_metadata[,2:11]

# Merge file info with metadata by data_id
all_metadata <- merge(clinical_metadata, all_metadata, by = "data_id")

# Trim filenames in metadata
all_metadata$filenames <- sapply(all_metadata$filenames, FUN = function(x) strsplit(x, "/")[[1]][[4]])

################################################################################
# Check for missing data
################################################################################

# Entries metadata CSV that don't have files
clinical_metadata$data_id[!clinical_metadata$data_id %in% all_metadata$data_id]

# Files that aren't in metadata csv
names(final_flowset@frames)[!names(final_flowset@frames) %in% all_metadata$filenames]

################################################################################
# Only retain files/samples that appear in both flowset and metadata
################################################################################

# Filter out non-blood samples
all_metadata <- all_metadata[all_metadata$sample_type != "pleural_effusion",]

files_to_keep <- intersect(all_metadata$filenames,names(final_flowset@frames))

final_flowset <- final_flowset[files_to_keep]


metadata_to_use <- all_metadata %>% 
  dplyr::filter(filenames %in% files_to_keep)

################################################################################
# Generate SingleCellExperiment object
################################################################################
#Get colnames for metadata factors
factor_colnames <- colnames(metadata_to_use)[-which(colnames(metadata_to_use) == "filenames")]
factor_colnames <- factor_colnames[-which(colnames(metadata_to_use) == "data_id")]

# construct SingleCellExperiment
sce <- prepData(final_flowset, panel = marker_info, md=metadata_to_use, features = marker_info$fcs_colname,
                md_cols = list(file = "filenames", id = "data_id", factors = factor_colnames),
                transform = T)

# sce <- readRDS("data/cytof_objects/sclc_all_samples_object.rds")

all_experiments <- colData(sce) %>% 
  as.data.frame() %>% 
  count(experiment_id, condition,sample_type) %>% 
  count(experiment_id) %>% 
  filter(n > 2) %>% 
  pull(experiment_id) %>% 
  as.character()


sce <- sce[,sce$experiment_id %in% all_experiments]

sce$cell_id <- paste0("cell_",1:ncol(sce))

ctc_ids <- c()


curr_experiment <- all_experiments[1]
for(curr_experiment in all_experiments){
  curr_sce <- sce[,sce$experiment_id == curr_experiment]
  unique(curr_sce$sample_type)
  
  is_cancer_patient <- colData(curr_sce)$condition == "cancer"
  is_normal_patient <- colData(curr_sce)$condition == "normal"
  
  dna_channel <- "Ir191"
  
  dna_values <- assay(curr_sce, "exprs")[dna_channel, ]
  dna_normal <- dna_values[is_normal_patient]
  
  # Robust stats: use median and MAD
  median_dna <- median(dna_normal, na.rm = TRUE)
  mad_dna <- mad(dna_normal, na.rm = TRUE)
  
  # Define a normal range (e.g., ±3 MADs)
  lower_thresh <- median_dna - 2 * mad_dna
  upper_thresh <- median_dna + 2 * mad_dna
  
  
  is_aneuploid <- dna_values < lower_thresh | dna_values > upper_thresh
  colData(curr_sce)$aneuploid <- is_aneuploid
  
  curr_ctcs <- curr_sce$cell_id[colData(curr_sce)$aneuploid == T &
                                  colData(curr_sce)$condition == "cancer" &
                                  colData(curr_sce)$sample_type == "blood"]
  
  ctc_ids <- append(ctc_ids,curr_ctcs)
  
}



length(ctc_ids)

(length(ctc_ids)/ncol(sce))*100



ctcs <- sce[,sce$cell_id %in% ctc_ids]





