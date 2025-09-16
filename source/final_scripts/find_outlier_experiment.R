source("source/sclc_cytof_functions.R")

set.seed(42)
################################################################################
# Read in data
################################################################################
sce <- readRDS("data/cytof_objects/sclc_all_samples_object_no_qc.rds")

# Get samples that are run in multiple experiments
to_test <- as.data.frame(sce@colData) %>%
  dplyr::select(experiment_id,collection_id) %>%
  distinct() %>%
  dplyr::count(collection_id) %>%
  arrange(desc(n)) %>%
  dplyr::filter(n > 1) %>%
  pull(collection_id) %>%
  as.character()

samples_with_outlier_exp <- as.data.frame(sce@colData) %>%
  dplyr::select(experiment_id,collection_id) %>% 
  filter(experiment_id== "531050") %>% 
  pull(collection_id) %>% 
  unique() %>%
  as.character()

intersect(samples_with_outlier_exp,to_test)



sort(to_test)

# Checking SC454-1
temp <- sce[markers,sce$collection_id == "SC454-1"]
plotExprs(temp, color_by = "experiment_id", assay = "exprs")

temp <- sce[markers,sce$collection_id == "SC435-1"]
plotExprs(temp, color_by = "experiment_id", assay = "exprs")

temp <- sce[markers,sce$collection_id == "SC442-2"]
plotExprs(temp, color_by = "experiment_id", assay = "exprs")

temp <- sce[markers,sce$collection_id == "NORMAL8-1"]
plotExprs(temp, color_by = "experiment_id", assay = "exprs")

temp <- sce[markers,sce$collection_id == "NORMAL7-1"]
plotExprs(temp, color_by = "experiment_id", assay = "exprs")

temp <- sce[markers,sce$collection_id == "SC442-3"]
plotExprs(temp, color_by = "experiment_id", assay = "exprs")

temp <- sce[markers,sce$collection_id == "SC454-7"]
plotExprs(temp, color_by = "experiment_id", assay = "exprs")




