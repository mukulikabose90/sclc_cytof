
sce <- readRDS("data/cytof_objects/sclc_all_samples_object.rds")

sample_df <- list()

################################################################################
# number of total liquid biopsies
sce %>% 
  colData() %>%
  as.data.frame() %>% 
  dplyr::filter(condition == "cancer") %>% 
  pull(collection_id) %>% 
  unique() %>% 
  as.character() %>% 
  length()

################################################################################
# normal patients
count1 <- sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  dplyr::filter(condition == "normal") %>% 
  pull(patient_id) %>% 
  unique() %>% 
  gsub("NORMAL","",.) %>% 
  as.numeric() %>% 
  sort() %>% 
  paste0("NORMAL",.) %>% 
  length()

# normal LBs
count2 <- sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  dplyr::filter(condition == "normal") %>% 
  pull(collection_id) %>% 
  unique() %>% 
  length()

sample_df[["Patients"]] <- append(sample_df[["Patients"]], count1)
sample_df[["Liquid Biopsies"]] <- append(sample_df[["Liquid Biopsies"]], count2)

################################################################################
# cancer patients
count1 <- sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  dplyr::filter(condition == "cancer") %>% 
  pull(patient_id) %>% 
  unique() %>% 
  gsub("SC","",.) %>% 
  as.numeric() %>% 
  sort() %>% 
  paste0("SC",.) %>% 
  length()


# cancer patient LBs
count2 <- sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  dplyr::filter(condition == "cancer") %>% 
  pull(collection_id) %>% 
  unique() %>% 
  length()

sample_df[["Patients"]] <- append(sample_df[["Patients"]], count1)
sample_df[["Liquid Biopsies"]] <- append(sample_df[["Liquid Biopsies"]], count2)
################################################################################
# naive patients
count1 <- sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  dplyr::filter(treatment_status == "naive") %>% 
  pull(patient_id) %>% 
  unique() %>% 
  length()

# naive patient LBs
count2 <- sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  dplyr::filter(treatment_status == "naive") %>% 
  pull(collection_id) %>% 
  unique() %>% 
  length()

sample_df[["Patients"]] <- append(sample_df[["Patients"]], count1)
sample_df[["Liquid Biopsies"]] <- append(sample_df[["Liquid Biopsies"]], count2)

################################################################################
# treated with SOC patients
count1 <- sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  dplyr::filter(treatment_status == "treated" & (tarla == "pre" | is.na(tarla))) %>% 
  pull(patient_id) %>% 
  unique() %>% 
  length()


# treated with SOC patient LBs
count2 <- sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  dplyr::filter(treatment_status == "treated" & (tarla == "pre" | is.na(tarla))) %>% 
  pull(collection_id) %>% 
  unique() %>% 
  length()

sample_df[["Patients"]] <- append(sample_df[["Patients"]], count1)
sample_df[["Liquid Biopsies"]] <- append(sample_df[["Liquid Biopsies"]], count2)

################################################################################
# longitudinal getting SOC patients
soc_long_pts <- sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  dplyr::filter(condition == "cancer" & is.na(tarla)) %>% 
  count(patient_id,collection_id) %>% 
  count(patient_id) %>% 
  filter(n > 1) %>% 
  pull(patient_id) %>% 
  as.character()

count1 <- length(soc_long_pts)

# longitudinal getting SOC LBs
count2 <- sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  dplyr::filter(patient_id %in% soc_long_pts) %>% 
  pull(collection_id) %>% 
  unique() %>% 
  length()

sample_df[["Patients"]] <- append(sample_df[["Patients"]], count1)
sample_df[["Liquid Biopsies"]] <- append(sample_df[["Liquid Biopsies"]], count2)
  
################################################################################
#  tarla patients
count1 <- sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  dplyr::filter(tarla == "pre" | tarla == "post") %>% 
  pull(patient_id) %>% 
  unique() %>% 
  as.character() %>% 
  length()

# pre tarla patient LBs
count2 <- sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  dplyr::filter(tarla == "pre" | tarla == "post") %>% 
  pull(collection_id) %>% 
  unique() %>% 
  as.character() %>% 
  length()

sample_df[["Patients"]] <- append(sample_df[["Patients"]], count1)
sample_df[["Liquid Biopsies"]] <- append(sample_df[["Liquid Biopsies"]], count2)

################################################################################
# pre tarla patients
count1 <- sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  dplyr::filter(tarla == "pre") %>% 
  pull(patient_id) %>% 
  unique() %>% 
  as.character() %>% 
  length()

# pre tarla patient LBs
count2 <- sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  dplyr::filter(tarla == "pre") %>% 
  pull(collection_id) %>% 
  unique() %>% 
  as.character() %>% 
  length()

sample_df[["Patients"]] <- append(sample_df[["Patients"]], count1)
sample_df[["Liquid Biopsies"]] <- append(sample_df[["Liquid Biopsies"]], count2)

################################################################################
# post tarla patients
count1 <- sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  dplyr::filter(tarla == "post") %>% 
  pull(patient_id) %>% 
  unique() %>% 
  as.character() %>% 
  length()

# post tarla patient LBs
count2 <- sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  dplyr::filter(tarla == "post") %>% 
  pull(collection_id) %>% 
  unique() %>% 
  as.character() %>% 
  length()

sample_df[["Patients"]] <- append(sample_df[["Patients"]], count1)
sample_df[["Liquid Biopsies"]] <- append(sample_df[["Liquid Biopsies"]], count2)

################################################################################
# longitudinal getting tarla patients
tarla_long_pts <- sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  dplyr::filter(condition == "cancer" & !is.na(tarla)) %>% 
  count(patient_id,collection_id) %>% 
  count(patient_id) %>% 
  filter(n > 1) %>% 
  pull(patient_id) %>% 
  as.character()

count1 <- length(soc_long_pts)

# longitudinal getting tarla LBs
count2 <- sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  dplyr::filter(patient_id %in% tarla_long_pts) %>% 
  pull(collection_id) %>% 
  unique() %>% 
  length()

sample_df[["Patients"]] <- append(sample_df[["Patients"]], count1)
sample_df[["Liquid Biopsies"]] <- append(sample_df[["Liquid Biopsies"]], count2)

################################################################################
final_df <- as.data.frame(sample_df)

rownames(final_df) <- c("Healthy Donors","SCLC Samples","Naive","Treated SOC","Longitudinal SOC","Tarla","Pre-Tarla","Post-Tarla","Longitudinal Tarla")

final_df





intersect(pre_tarla_pt,post_tarla_pt)
