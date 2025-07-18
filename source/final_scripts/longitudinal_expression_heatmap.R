source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################

markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam", "E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist", "SLUG", "PD-L1", "p-YAP", "CD44", "CD24")

col_fun = colorRamp2(c(-2, -1, 0, 1, 2), 
                     c("#313695",  # deep blue
                       "#74add1",  # light blue
                       "#f7f7f7",  # white (center, 0)
                       "#f46d43",  # light red
                       "#a50026"))

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

long_patients <- as.data.frame(ctcs@colData) %>% 
  filter(tarla != "post" | is.na(tarla)) %>% 
  select(patient_id,sample_num) %>% 
  distinct() %>% 
  dplyr::count(patient_id) %>% 
  dplyr::filter(n > 1) %>% 
  pull(patient_id) %>% 
  as.character()

patients_to_remove <- ctcs %>% 
  colData() %>% 
  as.data.frame() %>% 
  count(patient_id,sample_num) %>% 
  filter(n < 10) %>% 
  pull(patient_id) %>% 
  unique() %>% 
  as.character()

long_patients <- long_patients[!long_patients %in% patients_to_remove]


long_data <- ctcs[,ctcs$patient_id %in% long_patients]

#Scale expression
assay(long_data, "exprs") <- t(scale(t(assay(long_data, "exprs"))))

heatmap_list <- list()
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
  curr_ht <- Heatmap(curr_heatmap, cluster_rows = F, cluster_columns = F,column_names_rot = 0,
                     name="Median Scaled\n   Expression", column_title = curr_patient, col = col_fun)
  
  
  heatmap_list <- append(heatmap_list, list(curr_ht))
  
}



dummy_ht <- draw(heatmap_list[[1]])

legend_obj <- color_mapping_legend(dummy_ht@ht_list[[1]]@matrix_color_mapping, plot = FALSE)

jpeg(glue("figures/longitudinal_expression_heatmap.jpg"), width=260,height=150, units = "mm", res=1000)

# Now draw the 16 heatmaps without legend in a 4x4 grid
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 6, widths = unit.c(unit(rep(1, 5), "null"), unit(2.8, "cm")))))


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

# Add the legend on the right side of the page
pushViewport(viewport(x = 0.95, y = 0.5, width = 1, height = 1, just = c("right", "center"), layout.pos.col = 6))
grid.draw(legend_obj)
upViewport()


dev.off()

################################################################################
# Remove tarla patients and plot heatmaps again
# patients_to_use <- ctcs %>% 
#   colData() %>% 
#   as.data.frame() %>% 
#   filter(tarla != "post") %>% 
#   select(patient_id,collection_id) %>% 
#   distinct() %>% 
#   count(patient_id) %>% 
#   filter(n > 1) %>% 
#   pull(patient_id)
# 
# 
# 
# length(heatmap_list)
# 
# dummy_ht <- draw(heatmap_list[[1]])
# 
# legend_obj <- color_mapping_legend(dummy_ht@ht_list[[1]]@matrix_color_mapping, plot = FALSE)
# 
# jpeg(glue("figures/longitudinal_notarla_expression_heatmap.jpg"), width=300,height=100, units = "mm", res=1000)
# 
# # Now draw the 16 heatmaps without legend in a 4x4 grid
# grid.newpage()
# pushViewport(viewport(layout = grid.layout(1,6, widths = unit.c(unit(rep(1, 3), "null"), unit(2.8, "cm")))))
# pushViewport(viewport(layout = grid.layout(1,6, widths = unit.c(unit(rep(1, 5), "null"), unit(2.8, "cm")))))
# 
# for (i in seq_along(heatmap_list)) {
#   col <- i
#   vp <- viewport(layout.pos.row = 1, layout.pos.col = col)
#   pushViewport(vp)
#   draw(heatmap_list[[i]], newpage = FALSE,
#        show_heatmap_legend = FALSE,
#        show_annotation_legend = FALSE)
#   upViewport()
# }
# # for (i in seq_along(heatmap_list)) {
# #   row <- ceiling(i / 3)
# #   col <- (i - 1) %% 3 + 1
# #   vp <- viewport(layout.pos.row = row, layout.pos.col = col)
# #   pushViewport(vp)
# #   draw(heatmap_list[[i]], newpage = FALSE,
# #        show_heatmap_legend = FALSE,
# #        show_annotation_legend = FALSE)
# #   upViewport()
# # }
# 
# # Add the legend on the right side of the page
# pushViewport(viewport(x = 0.95, y = 0.5, width = 1, height = 1, just = c("right", "center"), layout.pos.col = 6))
# grid.draw(legend_obj)
# upViewport()
# 
# 
# dev.off()


