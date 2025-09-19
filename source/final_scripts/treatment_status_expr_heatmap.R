################################################################################
# This script plots a heatmap of the expression of all protein markers between
# naive CTCs and CTCs treated with SOC in patients that have > 10 cells in both
# treatment status
################################################################################
source("source/sclc_cytof_functions.R")

set.seed(42)
################################################################################
# Read in data
################################################################################
ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

################################################################################
# Set up colors and proteins
################################################################################
col_fun = colorRamp2(c(-2, -1, 0, 1, 2), 
                     c("#313695",  # deep blue
                       "#74add1",  # light blue
                       "#f7f7f7",  # white (center, 0)
                       "#f46d43",  # light red
                       "#a50026"))

markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam", "E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist", "SLUG", "PD-L1", "p-YAP", "CD44", "CD24")

################################################################################
# Select patients to use
################################################################################
patients_to_use <- ctcs %>% 
  colData() %>% 
  as.data.frame() %>% 
  filter(tarla != "post" | is.na(tarla)) %>% 
  select(patient_id,treatment_status) %>% 
  distinct() %>% 
  count(patient_id) %>% 
  filter(n > 1) %>% 
  pull(patient_id) %>% 
  as.character()

# Remove samples with too few cells
patients_to_remove <- ctcs %>% 
  colData() %>% 
  as.data.frame() %>% 
  count(patient_id,sample_num) %>% 
  filter(n < 10) %>% 
  pull(patient_id) %>% 
  unique() %>% 
  as.character()

patients_to_use <- patients_to_use[!patients_to_use %in% patients_to_remove]

curr_data <- ctcs[,ctcs$patient_id %in% patients_to_use]

################################################################################
# Plot marker expression between naive and treated (for each patient)
################################################################################

#Scale expression
assay(curr_data, "exprs") <- t(scale(t(assay(curr_data, "exprs"))))

all_ht <- list()
for(curr_patient in patients_to_use){
  
  curr_heatmap <- list()
  
  for(i in c("naive","treated")){
    
    y <- assay(curr_data[,curr_data$patient_id == curr_patient & curr_data$treatment_status == i], "exprs")
    
    y <- y[markers_to_use,]
    
    curr_heatmap[[i]] <- rowMedians(y)
    
  }
  
  curr_heatmap <- do.call(cbind, curr_heatmap)
  
  colnames(curr_heatmap) <- c("Naive","Treated")
  rownames(curr_heatmap) <- markers_to_use
  curr_ht <- Heatmap(curr_heatmap, cluster_rows = F, cluster_columns = F,column_names_rot = 45,
                     name="Median Scaled\n   Expression", column_title = curr_patient, col = col_fun)
  
  all_ht <- append(all_ht, list(curr_ht))
}

################################################################################
# Save figure
################################################################################
tiff(glue("figures/treatment_status_expression_heatmap.tiff"), width=160,height=160, units = "mm", res=600)
print(all_ht[[1]])
dev.off()




