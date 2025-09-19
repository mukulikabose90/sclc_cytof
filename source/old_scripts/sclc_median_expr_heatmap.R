source("source/cytof_de_function.R")

scale_values <- function(x){(x-min(x))/(max(x)-min(x))}

col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
################################################################################

sce <- readRDS("data/cytof_objects/sclc_cytof_sce_object.rds")

sce <- sce[,colData(sce)$sample_type != "cell_line"]

sce <- sce[rowData(sce)$marker_class == "state",]

experiments <- levels(colData(sce)$experiment_id)
# 
# curr_experiment <- experiments[2]
# 
# curr_exp_sce <- sce[,colData(sce)$experiment_id == curr_experiment]
# 
curr_exp_sce <- sce[,colData(sce)$experiment_id %in% experiments]


y <- assay(curr_exp_sce, "exprs")

df <- data.frame(t(y), colData(curr_exp_sce), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")
gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(curr_exp_sce)))

################################################################################
# Patient Heatmap
################################################################################
temp <- gg_df %>% 
  group_by(patient_id, antigen) %>% 
  summarise(mean(expression)) %>% 
  pivot_wider(names_from = "antigen", values_from = "mean(expression)") %>% 
  column_to_rownames("patient_id")

condition_colors <- as.data.frame(colData(curr_exp_sce)) %>% 
  select(patient_id,condition) %>% 
  distinct() %>% 
  deframe() 

patient_names <- names(condition_colors)
names(condition_colors) <- patient_names

ha = rowAnnotation("Sample Type" = ifelse(grepl("SC", rownames(temp)), "Cancer","Normal"), 
                   col = list("Sample Type"=c("Cancer" = "lightgreen","Normal" = "pink")),
                   show_annotation_name = F)

temp_scaled <- apply(temp, MARGIN = 2, FUN = function(x) scale_values(x))


# Heatmap(temp_scaled, name = "Mean\nExpression", right_annotation = ha, col = col_fun)

################################################################################
# Sample Heatmap
################################################################################
temp <- gg_df %>% 
  group_by(id, antigen) %>% 
  summarise(mean(expression)) %>% 
  pivot_wider(names_from = "antigen", values_from = "mean(expression)") %>% 
  column_to_rownames("id")

condition_colors <- as.data.frame(colData(sce)) %>% 
  select(id,condition) %>% 
  distinct() %>% 
  deframe() 

patient_names <- names(condition_colors)
names(condition_colors) <- patient_names


heatmap_metadata <- data.frame(colData(sce)) %>% 
  dplyr::select(id,patient_id,condition,experiment_id) %>% 
  distinct() 



all(heatmap_metadata$sample_id == rownames(temp))

# Create Annotations
condition_anno <-  rowAnnotation("Condition" = as.vector(heatmap_metadata$condition), 
                   col = list("Condition"=c("cancer" = "lightgreen","normal" = "pink")),
                   show_annotation_name = T)

ttxx <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(heatmap_metadata$patient_id)))
names(ttxx) <- unique(heatmap_metadata$patient_id)
patient_anno <- rowAnnotation("Patient ID" = heatmap_metadata$patient_id, 
                   col = list("Patient ID"= ttxx),
                   show_annotation_name = T)


ttxx <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(heatmap_metadata$experiment_id)))
names(ttxx) <- unique(heatmap_metadata$experiment_id)
experiment_anno <- rowAnnotation("Experiment" = heatmap_metadata$experiment_id, 
                             col = list("Experiment"= ttxx),
                             show_annotation_name = T)


dim(temp)

# temp_scaled <- apply(temp, MARGIN = 2, FUN = function(x) scale_values(x))
temp_scaled <- scale(temp)



Heatmap(temp_scaled, name = "Mean\nExpression", right_annotation = c(condition_anno, patient_anno, experiment_anno),
        cluster_rows = T, cluster_columns = T)



