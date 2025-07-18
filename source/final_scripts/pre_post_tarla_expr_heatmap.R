source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)

col_fun = colorRamp2(c(-2, -1, 0, 1, 2), 
                     c("#313695",  # deep blue
                       "#74add1",  # light blue
                       "#f7f7f7",  # white (center, 0)
                       "#f46d43",  # light red
                       "#a50026"))

################################################################################

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")


patients_to_use <- ctcs %>% 
  colData() %>% 
  as.data.frame() %>% 
  filter(!is.na(tarla)) %>% 
  select(patient_id,tarla) %>% 
  distinct() %>% 
  count(patient_id) %>% 
  filter(n > 1) %>% 
  pull(patient_id)


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


tarla_ctcs <- ctcs[,ctcs$patient_id %in% patients_to_use]

#Scale expression
assay(tarla_ctcs, "exprs") <- t(scale(t(assay(tarla_ctcs, "exprs"))))


all_ht <- list()
for(curr_patient in patients_to_use){
  
  
  curr_heatmap <- list()
  
  for(i in c("pre","post")){
    
    y <- assay(tarla_ctcs[,tarla_ctcs$patient_id == curr_patient & tarla_ctcs$tarla == i], "exprs")
    
    y <- y[markers_to_use,]
    
    curr_heatmap[[i]] <- rowMedians(y)
    
  }
  
  curr_heatmap <- do.call(cbind, curr_heatmap)
  
  colnames(curr_heatmap) <- c("Pre-Tarlatamab","Post-Tarlatamab")
  rownames(curr_heatmap) <- markers_to_use
  curr_ht <- Heatmap(curr_heatmap, cluster_rows = F, cluster_columns = F,column_names_rot = 45,
          name="Median Scaled\n   Expression", column_title = curr_patient, col = col_fun)
  
  all_ht <- append(all_ht, list(curr_ht))
}


jpeg(glue("figures/pre_post_tarla_expression_heatmap.jpg"), width=300,height=160, units = "mm", res=1000)
print(all_ht[[1]]+all_ht[[2]]+all_ht[[3]])
dev.off()


