# This script calculates the median expression for each protein for each
# patient. It then plots a scaled heatmap.

################################################################################
library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(CATALYST)
library(reshape2)
source("source/cytof_de_function.R")

col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))



set.seed(42)
################################################################################
# Read in CyTOF data with cluster assignments
sce <- readRDS("data/cytof_objects/sclc_cytof_sce_object_with_subtypes.rds")

curr_exp <- levels(colData(sce)$experiment_id)[3]

# sce <- sce[,colData(sce)$experiment_id == curr_exp]


sce <- sce[,colData(sce)$condition == "cancer"]

# Subset to cell in cluster 1
ctcs <- sce[,colData(sce)$new_clusters == 1]

# Subset to only state markers
ctcs <- ctcs[rowData(ctcs)$marker_class == "state",]

# Select expression data
y <- assay(ctcs, "exprs")

#Create tidy dataframe for each sample
df <- data.frame(t(y), colData(ctcs), check.names = FALSE)
value <- ifelse("exprs" == "exprs", "expression", "exprs")
gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(ctcs)))

# Create heatmap of median expression for each patient
heatmap <- gg_df %>% 
  group_by(patient_id, sample_num, antigen) %>% 
  summarise(median(expression)) %>% 
  pivot_wider(names_from = "antigen", values_from = "median(expression)") %>% 
  mutate(sample_num_id = paste0(patient_id,"-",sample_num)) %>% 
  column_to_rownames("sample_num_id") %>% 
  select(-c(patient_id,sample_num))

heatmap_scaled <-  apply(heatmap, MARGIN = 2, FUN = function(x) scale_values(x))

##################################################################

heatmap_metadata <- data.frame(colData(ctcs)) %>% 
  dplyr::select(sample_id,patient_id,condition,sample_num, treatment_status,subtypes) %>% 
  mutate(sample_num_id = paste0(patient_id,"-",sample_num)) %>% 
  select(-c(sample_id)) %>% 
  distinct()

all(heatmap_metadata$sample_id == rownames(heatmap_scaled))

# Create Annotations

ttxx <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(heatmap_metadata$patient_id)))
names(ttxx) <- unique(heatmap_metadata$patient_id)
patient_anno <- rowAnnotation("Patient ID" = heatmap_metadata$patient_id, 
                              col = list("Patient ID"= ttxx),
                              show_annotation_name = T)


ttxx <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(heatmap_metadata$sample_num)))
names(ttxx) <- unique(heatmap_metadata$sample_num)
sample_num_anno <- rowAnnotation("Sample Number" = heatmap_metadata$sample_num, 
                                 col = list("Sample Number"= ttxx),
                                 show_annotation_name = T)

treatment_anno <-  rowAnnotation("Status" = as.vector(heatmap_metadata$treatment_status), 
                                 col = list("Status"=c("treated" = "darkorange","naive" = "royalblue")),
                                 show_annotation_name = T)

ttxx <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(heatmap_metadata$subtypes)))
names(ttxx) <- unique(heatmap_metadata$subtypes)
subtype_anno <- rowAnnotation("Subtype" = heatmap_metadata$subtypes, 
                                 col = list("Subtype"= ttxx),
                                 show_annotation_name = T)

Heatmap(heatmap_scaled, name = "Median\nExpression", right_annotation = c(subtype_anno,treatment_anno),
        cluster_rows = T, cluster_columns = T, col=col_fun)







