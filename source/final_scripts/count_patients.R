

sce <- readRDS("data/cytof_objects/sclc_all_samples_object.rds")



normal_patients <- sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  filter(condition == "normal") %>% 
  pull(patient_id) %>% 
  unique() %>% 
  gsub("NORMAL","",.) %>% 
  as.numeric() %>% 
  sort() %>% 
  paste0("NORMAL",.)

write.table(cancer_patients,"data/test.txt",quote = F,row.names = F)

cancer_patients <- sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  filter(condition == "cancer") %>% 
  pull(patient_id) %>% 
  unique() %>% 
  gsub("SC","",.) %>% 
  as.numeric() %>% 
  sort() %>% 
  paste0("SC",.)
