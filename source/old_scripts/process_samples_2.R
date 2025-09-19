source("source/cytof_de_function.R")

create_metadata <- function(files,experiment_id){
  
  # Meta-data: experiment information
  # experiment_id <- strsplit(experiment_id,"_")[[1]][2]
  
  filenames <- as.character(files)
  
  # sample information
  data_id <- gsub("_norm_CD45-CD56\\+", "", gsub("\\.fcs$", "", filenames))

  data_id <- sapply(data_id, FUN = function(x) strsplit(x,"/")[[1]][4])
  
  data_id <- sapply(data_id, FUN = function(x) strsplit(x,"_")[[1]][1])
  
  data_id <- paste0(data_id,"_",experiment_id)
  
  # experiment_metadata <- data.frame(cbind(filenames,data_id,sample_id,patient_id,sample_num,experiment_id))
  experiment_metadata <- data.frame(cbind(filenames,data_id))
  rownames(experiment_metadata) <- NULL
  
  return(experiment_metadata)
}

###################################################################################

all_experiments <- list.files("data/cytof_raw_data")

# original_experiments <- c("experiment_513549","experiment_514521","experiment_515600")
# new_experiments <- c("experiment_508095","experiment_508814","experiment_511467")
# all_experiments <- new_experiments

markers_to_use <- readRDS("data/markers_to_use.rds")

all_metadata <- list()
all_files <- c()
all_flowsets <- list()
for(curr_experiment in all_experiments){
  curr_experiment_id <- strsplit(curr_experiment, "_")[[1]][2]
  
  # Get fcs file names of current experiment
  curr_files <- list.files(glue("data/cytof_raw_data/{curr_experiment}"), pattern = "\\.fcs$", full.names = T)

  curr_flowSet <- read.flowSet(curr_files, transformation = FALSE, truncate_max_range = FALSE)
  
  # Subset flowSet to include only selected markers
  # Apply subsetting to each flowFrame
  curr_flowSet_subset <- fsApply(curr_flowSet, function(ff) {
    ff[, markers_to_use]
  })
  
  all_flowsets <- append(all_flowsets,curr_flowSet_subset)
  
  
  all_files <- append(all_files, curr_files)

  # Create metadata
  all_metadata <- append(all_metadata, list(create_metadata(curr_files,curr_experiment_id)))
}

final_flowset <- Reduce(rbind2, all_flowsets)


# add clinical metadata
all_metadata <- do.call(rbind, all_metadata)
clinical_metadata <- read.csv("data/cytof_metadata.csv")

# Select the clinical metadata you want to keep
clinical_metadata <- clinical_metadata[,1:12]

# all_metadata <- merge(clinical_metadata, all_metadata, by = colnames(all_metadata)[-1])
all_metadata <- merge(clinical_metadata, all_metadata, by = "data_id")


# Trim filnames in metadata
all_metadata$filenames <- sapply(all_metadata$filenames, FUN = function(x) strsplit(x, "/")[[1]][[4]])

# Subset flowset and metadata to only include blood samples
files_to_keep <- all_metadata %>% 
  dplyr::filter(sample_type == "blood") %>% 
  pull(filenames)

final_flowset <- final_flowset[files_to_keep]

metadata_to_use <- all_metadata%>% 
  dplyr::filter(filenames %in% files_to_keep)

################################################################################
# Meta-data: marker information

marker_info <- read.csv("data/cytof_panel_info.csv")

marker_info <- data.frame(marker_info, stringsAsFactors = FALSE)

colnames(marker_info)

################################################################################

#Get colnames for metadata factors
factor_colnames <- colnames(metadata_to_use)[-which(colnames(metadata_to_use) == "filenames")]
factor_colnames <- factor_colnames[-which(colnames(metadata_to_use) == "data_id")]

# construct SingleCellExperiment
sce <- prepData(final_flowset, panel = marker_info, md=metadata_to_use, features = marker_info$fcs_colname,
                md_cols = list(file = "filenames", id = "data_id", factors = factor_colnames))




batch <- as.factor(colData(sce)$experiment_id)

corrected_exprs <- removeBatchEffect(assay(sce, "exprs"), batch=batch)
assay(sce, "batch_corrected_exprs") <- corrected_exprs

corrected_counts <- removeBatchEffect(assay(sce, "counts"), batch=batch)
assay(sce, "batch_corrected_counts") <- corrected_counts

saveRDS(sce, "data/cytof_objects/sclc_cytof_sce_object.rds")

