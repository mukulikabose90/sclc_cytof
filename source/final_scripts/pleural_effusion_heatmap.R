################################################################################
# This script reads in all CyTOF FCS files along with the metadata associated 
# with each experiment. This also include pleural effusion sample from
# MDA-SC443 for comparison with blood samples
################################################################################
source("source/sclc_cytof_functions.R")

set.seed(42)
################################################################################

markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam", "E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist", "SLUG", "PD-L1", "p-YAP", "CD44", "CD24")

col_fun = colorRamp2(c(-2, -1, 0, 1, 2), 
                     c("#313695",  # deep blue
                       "#74add1",  # light blue
                       "#f7f7f7",  # white (center, 0)
                       "#f46d43",  # light red
                       "#a50026"))

create_metadata <- function(files,experiment_id){
  
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

curr_experiment <- all_experiments[1]
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
clinical_metadata <- clinical_metadata[,2:13]

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
pe_samples_to_keep <- c("SC443-1P")

# Filter out non-blood samples
all_metadata <- all_metadata[all_metadata$sample_type == "blood" | all_metadata$collection_id %in% pe_samples_to_keep,]
# all_metadata <- all_metadata[all_metadata$sample_type == "blood",]

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

sce$patient_id <- paste0("MDA-",sce$patient_id)
sce$collection_id <- paste0("MDA-",sce$collection_id)



sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  filter(patient_id %in% c("MDA-SC443")) %>% 
  pull(collection_id) %>% 
  unique()

samples_to_use <- c("MDA-SC443-1","MDA-SC443-1P")

curr_data <- sce

#Scale expression
assay(curr_data, "exprs") <- t(scale(t(assay(curr_data, "exprs"))))

curr_data <- curr_data[,curr_data$collection_id %in% samples_to_use]



markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam", "E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist", "SLUG", "PD-L1", "p-YAP", "CD44", "CD24")

curr_heatmap <- matrix(NA,ncol=2,nrow=length(markers_to_use))
for(curr_sample in samples_to_use){
  
  i <- which(curr_sample == samples_to_use)
  
  y <- assay(curr_data[,curr_data$collection_id == curr_sample], "exprs")
  
  y <- y[markers_to_use,]
  
  curr_heatmap[,i] <- rowMedians(y)
  
  
}

colnames(curr_heatmap) <- samples_to_use
rownames(curr_heatmap) <- markers_to_use

curr_ht <- Heatmap(curr_heatmap, cluster_rows = F, cluster_columns = F,column_names_rot = 0,column_names_centered = T,
                   name="Median Scaled\n   Expression", column_title = "", col = col_fun)


p1 <- draw(curr_ht)
p1

jpeg(glue("figures/pleural_effusion_expression_heatmap.jpg"), width=160,height=100, units = "mm", res=600)
print(p1)
dev.off()
