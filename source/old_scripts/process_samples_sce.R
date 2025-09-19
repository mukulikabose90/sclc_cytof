# This script reads in all SCLC CyTOF experiments after 9/12/24. It then compiles
# the metadata for the samples and creates a SCE object to contain all the data.
# The SCE object is then saved as an RDS file

library(tidyverse)
library(glue)
library(reshape2)
library(readxl)
library(HDCytoData)
library(CATALYST)
library(ComplexHeatmap)

# Functions to create patient and sample ID
create_patient_id <- function(x){
  if(grepl("NB",x)){ #Normal sample
    sample_id <- strsplit(x,"NB")[[1]][2]
    sample_id <- ifelse(is.na(sample_id),0,sample_id)
    patient_id <- glue("Normal_{sample_id}")
    
  } else if(grepl("SC",x)){ # Cancer samples
    
    #Create patient id
    
    
    patient_id <- strsplit(x, "-")[[1]][1]
    
    
  } else { #Cell line sample
    
    patient_id <- x
    
  }
  
  return(patient_id)
}

create_sample_num <- function(x){
  if(grepl("NB",x)){ #Normal sample
    
    sample_num <- strsplit(x,"NB")[[1]][2]
    sample_num <- ifelse(is.na(sample_num),1,sample_num)
    
    
  } else if(grepl("SC",x)){ # Cancer samples
    
    #Create  sample num
    
    sample_num <- strsplit(x, "-")[[1]][2]
    
  } else { #Cell line sample
    
    sample_num <- 1
    
  }
  
  return(sample_num)
}

create_metadata <- function(cytof_data_obj,experiment_id){
  # Create metadata dataframe
  patient_names <- sapply(names(cytof_data_obj@frames), function(x)strsplit(x,"_")[[1]][1])
  sample_id <- sapply(names(cytof_data_obj@frames), function(x)strsplit(x,"\\.")[[1]][1])
  
  condition <- ifelse(grepl("NB",patient_names), "normal", ifelse(grepl("SC",patient_names), "cancer", "cell_line"))
  patient_id <- sapply(patient_names, function(x) create_patient_id(x))
  sample_num <- sapply(patient_names, function(x) create_sample_num(x))
  
  metadata <- data.frame(names(cytof_data@frames),condition,patient_id,sample_id,experiment_id,sample_num)
  colnames(metadata)[1] <- "file_name"
  # rownames(metadata) <- NULL
  
  metadata$condition <- factor(metadata$condition, levels = c("normal", "cancer","cell_line"))
  metadata$sample_num <- factor(metadata$sample_num, 
                               levels = unique(metadata$sample_num[order(metadata$sample_id)]))
  metadata$experiment_id <- factor(metadata$experiment_id, 
                               levels = unique(metadata$experiment_id[order(metadata$sample_id)]))
  
  return(metadata)
}

markers_of_interest <- c("CD45", "CD56", "NeuroD1", "PD-L1", "SLUG", "Twist", "DLL3", "c-Casp3", "p-YAP",
                         "Vimentin", "CD44", "E-Cad", "MUC-1", "p-Rb", "Alcam", "CD24", "POU2F3","EpCAM","ASCL1")

################################################################################

all_experiments <- list.files("data/cytof_raw_data")

all_data <- list()
all_metadata <- list()
for(curr_experiment in all_experiments){
  experiment_id <- strsplit(curr_experiment, "_")[[1]][2]
  
  # Get fcs file names of current experiment
  fcs_file_names <- list.files(glue("data/cytof_raw_data/{curr_experiment}/"))
  
  # Create flowSet with current data
  cytof_data <- read.flowSet(files = fcs_file_names, path = glue("data/cytof_raw_data/{curr_experiment}/"))
  
  all_data <- append(all_data, (cytof_data))
  
  # Create metadata
  all_metadata <- append(all_metadata, list(create_metadata(cytof_data, experiment_id)))
}

cytof_data <- rbind2(all_data[[1]],all_data[[2]])

cytof_data <- rbind2(cytof_data, all_data[[3]])

all_metadata <- do.call(rbind, all_metadata)
rownames(all_metadata) <- NULL

# Create panel dataframe
fcs_colname <- colnames(cytof_data)

antigen <- sapply(pData(parameters(cytof_data[[1]]))$desc,FUN = function(x) strsplit(x,"_")[[1]][2]) 

sclc_panel <- as.data.frame(cbind(fcs_colname,antigen))
rownames(sclc_panel) <- NULL

sclc_panel <- sclc_panel %>% 
  dplyr::filter(!is.na(antigen)) %>% 
  dplyr::filter(antigen != "Dead" & antigen != "length")

sclc_panel$marker_class <- "none"

# Check that all panel markers are present in data
all(sclc_panel$fcs_colname %in% colnames(cytof_data))


# construct SingleCellExperiment
sce <- prepData(cytof_data, panel = sclc_panel, md=all_metadata, features = sclc_panel$fcs_colname,
                md_cols = list(file = "file_name", id = "sample_id", factors = c("condition", "patient_id","experiment_id","sample_num")))

# Filter to only include protein markers
sce <- sce[rownames(sce) %in% markers_of_interest,]
colData(sce)$cell_id <- 1:ncol(sce)

saveRDS(sce, "data/cytof_objects/sclc_cytof_sce_object.rds")








