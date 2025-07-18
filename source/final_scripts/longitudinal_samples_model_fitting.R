source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################

ctcs <- readRDS("data/cytof_objects/all_samples_ctcs_with_subtype.rds")

################################################################################
# remove collections with < 30 cells
################################################################################
samples_to_remove <- as.data.frame(ctcs@colData) %>%
  dplyr::count(collection_id) %>%
  dplyr::filter(n<5) %>%
  pull(collection_id) %>%
  as.character()

ctcs <- ctcs[,!ctcs$collection_id %in% samples_to_remove]

long_patients <- as.data.frame(ctcs@colData) %>% 
  select(patient_id,sample_num) %>% 
  distinct() %>% 
  dplyr::count(patient_id) %>% 
  dplyr::filter(n > 1) %>% 
  pull(patient_id) %>% 
  as.character()


long_data <- ctcs[,ctcs$patient_id %in% long_patients]

plot_df <- as.data.frame(long_data@colData) %>% 
  select(patient_id,sample_num,subtype,treatment_status) %>% 
  dplyr::count(patient_id,sample_num,subtype,treatment_status) %>% 
  group_by(patient_id,sample_num) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = n/total) %>% 
  mutate(sample_id = paste0(patient_id,"-",sample_num))

df <- plot_df %>% 
  dplyr::filter(subtype == "A") %>% 
  select(freq,sample_num,n,patient_id,treatment_status)


df$sample_num <- as.numeric(df$sample_num)

model <- glm("freq ~ sample_num + patient_id + treatment_status", data = df)

model <- glm("freq ~ treatment_status", data = df)

summary(model)


