source("../sclc_cytof/source/cytof_de_function.R")

create_metadata <- function(files, experiment_id){
  
  # Meta-data: experiment information
  
  filenames <- as.character(files)
  # # sample information
  data_id <- gsub("_norm_CD45-CD56\\+", "", gsub("\\.fcs$", "", filenames))
  # 
  treatment_status <- sapply(data_id, FUN = function(x) strsplit(x,"_")[[1]][1])
  timepoint <- sapply(data_id, FUN = function(x) strsplit(x,"_")[[1]][2])
  
  data_id <- gsub("_2025-0227","",data_id)
  experiment_metadata <- data.frame(cbind(filenames,data_id,timepoint,treatment_status))
  rownames(experiment_metadata) <- NULL
  
  experiment_metadata$day <- ifelse(treatment_status == "withdrawal", as.numeric(timepoint)+7, timepoint)
  
  return(experiment_metadata)
}

###################################################################################

all_experiments <- list.files("../H196_gemcitabine_treatment/data/raw_cytof_data")

curr_experiment <- strsplit(all_experiments, "_")[[1]][2]


curr_files <- list.files(glue("../H196_gemcitabine_treatment/data/raw_cytof_data/experiment_{curr_experiment}/"), pattern = "\\.fcs$", full.names = F)

# Create metadata
all_metadata <- create_metadata(curr_files, curr_experiment)

# add clinical metadata
# all_metadata <- do.call(rbind, all_metadata)
# clinical_metadata <- read.csv("data/cytof_metadata.csv")

# Select the clinical metadata you want to keep
# clinical_metadata <- clinical_metadata[,1:12]
# all_metadata <- merge(clinical_metadata, all_metadata, by = colnames(all_metadata)[-1])
# all_metadata <- merge(clinical_metadata, all_metadata, by = "data_id")
all_files <- list.files(glue("../H196_gemcitabine_treatment/data/raw_cytof_data/experiment_{curr_experiment}/"), pattern = "\\.fcs$", full.names = T)
d_flowSet <- read.flowSet(all_files, transformation = FALSE, truncate_max_range = FALSE)

# Filter data
# remove cell line samples before processing
# to_keep <- all_metadata$sample_type != "cell_line"s

# all_metadata <- all_metadata[to_keep,]


d_flowSet <- d_flowSet[all_metadata$filenames]

########################
# Meta-data: marker information

marker_info <- read.csv("data/cytof_panel_info.csv")

marker_info <- data.frame(marker_info, stringsAsFactors = FALSE)

colnames(marker_info)

################################################################################

#Get colnames for metadata factors
factor_colnames <- colnames(all_metadata)[-which(colnames(all_metadata) == "filenames")]
# factor_colnames <- factor_colnames[-which(colnames(all_metadata) == "data_id")]

# construct SingleCellExperiment
sce <- prepData(d_flowSet, panel = marker_info, md=all_metadata, features = marker_info$channel_name,
                md_cols = list(file = "filenames", id = "data_id", factors = factor_colnames))

# sce <- prepData(d_flowSet, panel = marker_info, md=all_metadata, features = marker_info$channel_name,
#                 md_cols = list(file = "filenames", id = "sample_id", factors = c("condition", "patient_id","sample_type", "experiment_id","sample_num")))



saveRDS(sce, "data/cytof_objects/drug_treatment_experiment_sce_object.rds")

