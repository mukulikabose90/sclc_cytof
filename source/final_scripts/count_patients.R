

sce <- readRDS("data/cytof_objects/sclc_all_samples_object.rds")





# normal patients
normal_patients <- sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  dplyr::filter(condition == "normal") %>% 
  pull(patient_id) %>% 
  unique() %>% 
  gsub("NORMAL","",.) %>% 
  as.numeric() %>% 
  sort() %>% 
  paste0("NORMAL",.)


length(normal_patients)
# write.table(cancer_patients,"data/test.txt",quote = F,row.names = F)

# cancer patients
cancer_patients <- sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  dplyr::filter(condition == "cancer") %>% 
  pull(patient_id) %>% 
  unique() %>% 
  gsub("SC","",.) %>% 
  as.numeric() %>% 
  sort() %>% 
  paste0("SC",.)

length(cancer_patients)

# naives
sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  dplyr::filter(treatment_status == "naive") %>% 
  distinct() %>% 
  arrange(date_run) %>% 
  select(collection_id) %>% 
  distinct() %>% 
  pull(collection_id) %>% 
  unique() %>% 
  as.character() %>% 
  length()

# treated
sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  dplyr::filter(treatment_status == "treated") %>% 
  distinct() %>% 
  arrange(date_run) %>% 
  select(collection_id) %>% 
  distinct() %>% 
  pull(collection_id) %>% 
  unique() %>% 
  as.character() %>% 
  length()


# number of liquid biopsies
sce %>% 
  colData() %>%
  as.data.frame() %>% 
  dplyr::filter(condition == "cancer") %>% 
  pull(collection_id) %>% 
  unique() %>% 
  as.character() %>% 
  length()


# longitudinal with tarla
sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  dplyr::filter(condition == "cancer") %>% 
  count(patient_id,collection_id) %>% 
  count(patient_id) %>% 
  filter(n > 1) %>% 
  nrow()

# longitudinal without tarla
sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  dplyr::filter(condition == "cancer" & is.na(tarla)) %>% 
  count(patient_id,collection_id) %>% 
  count(patient_id) %>% 
  filter(n > 1) %>% 
  nrow()


# pre tarla
pre_tarla_samples <- sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  dplyr::filter(tarla == "pre" ) %>% 
  pull(patient_id) %>% 
  unique() %>% 
  as.character() 

# post tarla
post_tarla_samples <- sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  dplyr::filter(tarla == "post") %>% 
  pull(patient_id) %>% 
  unique() %>% 
  as.character() 

intersect(sort(pre_tarla_samples),
          sort(post_tarla_samples))


sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  dplyr::filter(!is.na(tarla)) %>% 
  pull(patient_id) %>% 
  unique() %>% 
  as.character() %>% 
  length()

# Tarla-only Longitudinal 
sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  dplyr::filter(condition == "cancer" & !is.na(tarla)) %>% 
  count(patient_id,collection_id) %>% 
  count(patient_id) %>% 
  filter(n > 1) %>% 
  nrow()


