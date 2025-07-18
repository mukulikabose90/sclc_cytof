source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

################################################################################
# remove collections with < 30 cells
################################################################################
# samples_to_remove <- as.data.frame(ctcs@colData) %>%
#   dplyr::count(collection_id) %>%
#   dplyr::filter(n<30) %>%
#   pull(collection_id) %>%
#   as.character()
# 
# ctcs <- ctcs[,!ctcs$collection_id %in% samples_to_remove]



ctcs@colData %>% 
  as.data.frame() %>% 
  filter(treatment_status == "treated") %>% 
  pull(collection_id) %>% 
  unique() %>% 
  as.character() %>% 
  length()

ctcs@colData %>% 
  as.data.frame() %>% 
  filter(treatment_status == "naive") %>% 
  pull(collection_id) %>% 
  unique() %>% 
  as.character() %>% 
  length()

ctcs@colData %>% 
  as.data.frame() %>% 
  pull(collection_id) %>% 
  unique() %>% 
  as.character() %>% 
  length()

