################################################################################
# This script plots a heatmap of the expression of all protein markers in 
# each sample from longitudinal patientssthat have > 10 cells in both
# treatment status
################################################################################
source("source/sclc_cytof_functions.R")

set.seed(42)
################################################################################
# Read in data
################################################################################
ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

################################################################################
# Set up plot data
################################################################################
markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam", "E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist", "SLUG", "PD-L1", "p-YAP", "CD44", "CD24")

col_fun = colorRamp2(c(-2, -1, 0, 1, 2), 
                     c("#313695",  # deep blue
                       "#74add1",  # light blue
                       "#f7f7f7",  # white (center, 0)
                       "#f46d43",  # light red
                       "#a50026"))


long_patients <- as.data.frame(ctcs@colData) %>% 
  select(patient_id,sample_num) %>% 
  distinct() %>% 
  dplyr::count(patient_id) %>% 
  dplyr::filter(n > 1) %>% 
  pull(patient_id) %>% 
  as.character()

samples_to_remove <- ctcs %>% 
  colData() %>% 
  as.data.frame() %>% 
  count(patient_id,sample_num) %>% 
  filter(n < 10) %>% 
  mutate(collection_id = paste0(patient_id,"-",sample_num)) %>% 
  pull(collection_id) %>% 
  unique() %>% 
  as.character()


long_data <- ctcs[,ctcs$patient_id %in% long_patients & !ctcs$collection_id %in% samples_to_remove]

long_patients <- long_data@colData %>% 
  as.data.frame() %>% 
  count(patient_id,sample_num) %>% 
  count(patient_id) %>% 
  filter(n > 1) %>% 
  pull(patient_id)

#Scale expression
assay(long_data, "exprs") <- t(scale(t(assay(long_data, "exprs"))))


################################################################################
# Create heatmap for each patient
################################################################################
heatmap_list <- list()

long_patients <- c("MDA-SC392","MDA-SC430","MDA-SC293","MDA-SC338","MDA-SC500","MDA-SC454","MDA-SC501","MDA-SC547","MDA-SC254")

for(curr_patient in long_patients){
  
  curr_sample_ids <- long_data %>%
    colData() %>% 
    as.data.frame() %>% 
    filter(patient_id == curr_patient) %>% 
    pull(sample_num) %>% 
    unique() 
  
  
  curr_heatmap <- list()
  
  for(i in curr_sample_ids){
    
    sum(long_data$patient_id == curr_patient & long_data$sample_num == i)
    y <- assay(long_data[,long_data$patient_id == curr_patient & long_data$sample_num == i], "exprs")
    
    y <- y[markers_to_use,]
    
    curr_heatmap[[i]] <- rowMedians(y)
    
  }
  
  curr_heatmap <- do.call(cbind, curr_heatmap)
  
  colnames(curr_heatmap) <- curr_sample_ids
  colnames(curr_heatmap) <- 1:ncol(curr_heatmap)
  rownames(curr_heatmap) <- markers_to_use
  curr_ht <- Heatmap(curr_heatmap, cluster_rows = F, cluster_columns = F,column_names_rot = 0,column_names_centered = T,
                     name="Median Scaled\n   Expression", column_title = curr_patient, col = col_fun)
  
  heatmap_list <- append(heatmap_list, list(curr_ht))
  
}

dummy_ht <- draw(heatmap_list[[1]])

legend_obj <- color_mapping_legend(dummy_ht@ht_list[[1]]@matrix_color_mapping, plot = FALSE)

################################################################################
# Save figure
################################################################################

tiff(glue("figures/longitudinal_expression_heatmap.tiff"), width=260,height=150, units = "mm", res=600)

# Draw the heatmaps without legend in a grid
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 5, widths = unit.c(unit(rep(1, 5), "null"))))) #, unit(3, "cm")

for (i in seq_along(heatmap_list)) {
  row <- ceiling(i / 5)
  col <- (i - 1) %% 5 + 1
  vp <- viewport(layout.pos.row = row, layout.pos.col = col)
  pushViewport(vp)
  
  
  draw(heatmap_list[[i]], newpage = FALSE, 
       show_heatmap_legend = FALSE, 
       show_annotation_legend = FALSE)
  
  upViewport()
}

# Add the legend 
pushViewport(viewport(x = 0.95, y = 0, width = 1, height = 1, just = c("right", "center"), layout.pos.col = 5, layout.pos.row = 2))
grid.draw(legend_obj)
upViewport()


dev.off()

