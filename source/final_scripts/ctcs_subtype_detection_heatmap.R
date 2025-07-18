# This script plots a heatmap of the SCLC subtype TFs expression. Then assigns each cell
# to a heatmap based on the expression of the TFs.

################################################################################
source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################
# Read in CyTOF data with cluster assignments
################################################################################
cancer_enriched <- readRDS("data/cytof_objects/cancer_enriched_with_clusters.rds")


#Subset to clusters with known CTC marker profiles
# ctcs <- cancer_enriched[,cancer_enriched$new_clusters %in% c(2,5,6,7,8)]
ctcs <- cancer_enriched[,cancer_enriched$new_clusters %in% c(2,4,5,6,7,8)]

################################################################################
# Create heatmap
################################################################################

# Select expression data
y <- assay(ctcs, "exprs")

#Create tidy dataframe for each sample
df <- data.frame(t(y), colData(ctcs), check.names = FALSE)
value <- ifelse("exprs" == "exprs", "expression", "exprs")
gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(ctcs)))

# Create heatmap of protein expression for each cell
heatmap <- gg_df %>% 
  select(cell_id,antigen,expression) %>% 
  pivot_wider(names_from = antigen, values_from = expression) %>% 
  column_to_rownames("cell_id")

# Standardize expression within each protein
heatmap_scaled <- scale(heatmap)

# Subset to only SCLC subtype TFs 
heatmap_tf_scaled <- heatmap_scaled[,colnames(heatmap_scaled) %in% c("POU2F3","NeuroD1","ASCL1")]

#Set up metadata table to create annotations
heatmap_metadata <- data.frame(colData(ctcs)) %>% 
  dplyr::select(cell_id,collection_id,sample_id,patient_id,condition,sample_num, treatment_status,experiment_id) %>% 
  select(-c(sample_id)) %>% 
  distinct()

all_samples_heatmap <- t(heatmap_tf_scaled)
naive_heatmap <- all_samples_heatmap[,colnames(all_samples_heatmap) %in% ctcs[,ctcs$treatment_status == "naive"]$cell_id]
treated_heatmap <- all_samples_heatmap[,colnames(all_samples_heatmap) %in% ctcs[,ctcs$treatment_status == "treated"]$cell_id]

#############################################################################
# Create heatmap with all samples
#############################################################################
# col_fun = colorRamp2(c(-2, 0, 2), c("royalblue4","lightgray", "firebrick4"))
col_fun = colorRamp2(c(-2, -1, 0, 1, 2), 
                     c("#313695",  # deep blue
                       "#74add1",  # light blue
                       "#f7f7f7",  # white (center, 0)
                       "#f46d43",  # light red
                       "#a50026"))
# Create experiment ID annotation
# colors_to_use <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(heatmap_metadata$experiment_id)))
# names(colors_to_use) <- unique(heatmap_metadata$experiment_id)
# experiment_anno <- HeatmapAnnotation("Experiment ID" = heatmap_metadata$experiment_id, 
#                                  col = list("Experiment ID"= colors_to_use),
#                                  show_annotation_name = T,
#                                  annotation_legend_param = list(ncol=3))

# Create collection ID annotation
colors_to_use <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(heatmap_metadata$collection_id)))
names(colors_to_use) <- unique(heatmap_metadata$collection_id)
sample_anno <- HeatmapAnnotation("Sample ID" = heatmap_metadata$collection_id, 
                                 col = list("Sample ID"= colors_to_use),
                                 show_annotation_name = T,
                                 annotation_legend_param = list(ncol=3))

# Create patient ID annotation
colors_to_use <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(heatmap_metadata$patient_id)))
names(colors_to_use) <- unique(heatmap_metadata$patient_id)
patient_anno <- HeatmapAnnotation("Patient ID" = heatmap_metadata$patient_id, 
                                  col = list("Patient ID"= colors_to_use),
                                  show_annotation_name = T,
                                  annotation_legend_param = list(ncol=2))

# Get optimal k for kmeans
# silh_res <- fviz_nbclust(t(all_samples_heatmap), kmeans, method='silhouette')
# num_subclusters <- which(silh_res$data$y == max(silh_res$data$y))
# num_subclusters
num_subclusters <- 4

# col_fun = colorRamp2(c(min(heatmap_scaled), -1, 0, 1, max(heatmap_scaled)), 
#                      c("#313695",  # deep blue
#                        "#74add1",  # light blue
#                        "#f7f7f7",  # white (center, 0)
#                        "#f46d43",  # light red
#                        "#a50026"))



ht <- Heatmap(all_samples_heatmap, column_km = num_subclusters, top_annotation = sample_anno, name="Expression",
              cluster_columns = T, cluster_rows = F, show_column_names=F,col = col_fun,
              row_dend_reorder = F, column_title = c("SCLC-A", "SCLC-P","SCLC-I","SCLC-N"))

all_samples_ht <- draw(ht)

#############################################################################
# Create heatmap with naive samples only
#############################################################################
curr_metadata <- heatmap_metadata[heatmap_metadata$cell_id %in% colnames(naive_heatmap),]

# Create collection ID annotation
colors_to_use <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(curr_metadata$collection_id)))
names(colors_to_use) <- unique(curr_metadata$collection_id)
sample_anno <- HeatmapAnnotation("Sample ID" = curr_metadata$collection_id, 
                                 col = list("Sample ID"= colors_to_use),
                                 show_annotation_name = T,
                                 annotation_legend_param = list(ncol=3))

# Create patient ID annotation
colors_to_use <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(curr_metadata$patient_id)))
names(colors_to_use) <- unique(curr_metadata$patient_id)
patient_anno <- HeatmapAnnotation("Patient ID" = curr_metadata$patient_id, 
                                  col = list("Patient ID"= colors_to_use),
                                  show_annotation_name = T,
                                  annotation_legend_param = list(ncol=2))


ht <- Heatmap(naive_heatmap, column_km = num_subclusters, top_annotation = sample_anno, name="Expression",
              cluster_columns = T, cluster_rows = F, show_column_names=F,col = col_fun,
              row_dend_reorder = F, column_title = c("SCLC-N", "SCLC-P","SCLC-I","SCLC-A"))

naive_ht <- draw(ht)

#############################################################################
# Create heatmap with treated samples only
#############################################################################
curr_metadata <- heatmap_metadata[heatmap_metadata$cell_id %in% colnames(treated_heatmap),]

# Create collection ID annotation
colors_to_use <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(curr_metadata$collection_id)))
names(colors_to_use) <- unique(curr_metadata$collection_id)
sample_anno <- HeatmapAnnotation("Sample ID" = curr_metadata$collection_id, 
                                 col = list("Sample ID"= colors_to_use),
                                 show_annotation_name = T,
                                 annotation_legend_param = list(ncol=3))

# Create patient ID annotation
colors_to_use <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(curr_metadata$patient_id)))
names(colors_to_use) <- unique(curr_metadata$patient_id)
patient_anno <- HeatmapAnnotation("Patient ID" = curr_metadata$patient_id, 
                                  col = list("Patient ID"= colors_to_use),
                                  show_annotation_name = T,
                                  annotation_legend_param = list(ncol=2))



ht <- Heatmap(treated_heatmap, column_km = num_subclusters, top_annotation = sample_anno, name="Expression",
              cluster_columns = T, cluster_rows = F, show_column_names=F,col = col_fun,
              row_dend_reorder = F, column_title = c("SCLC-P", "SCLC-A","SCLC-I","SCLC-N"))

treated_ht <- draw(ht)

#############################################################################
# Find optimal number of subclusters
#############################################################################
p1 <- fviz_nbclust(t(all_samples_heatmap), kmeans, method='silhouette')+
  ggtitle("Optimal Number of Subclusters (CTCs)")

jpeg("figures/all_samples_ctcs_optimal_clusters.jpg", width = 200, height = 100, units = "mm", res = 1200)
print(p1)
dev.off()

#############################################################################
# ALL SAMPLES
#############################################################################
jpeg("figures/all_samples_ctcs_subtype_heatmap.jpg", width=300,height=120, units = "mm", res=1000)
print(all_samples_ht)
dev.off()
#############################################################################
# NAIVE SAMPLES
#############################################################################
jpeg("figures/all_samples_ctcs_naive_subtype_heatmap.jpg", width=300,height=100, units = "mm", res=1000)
print(naive_ht)
dev.off()
#############################################################################
#TREATED SAMPLES
#############################################################################
jpeg("figures/all_samples_ctcs_treated_subtype_heatmap.jpg", width=300,height=100, units = "mm", res=1000)
print(treated_ht)
dev.off()

#############################################################################
# Assign subtype to each cell based on kmeans clustering
#############################################################################
# Get column row from heatmap
clusters <- column_order(all_samples_ht)

# Rename clusters
names(clusters) <- c("A","P","I","N")
# names(clusters) <- c("P","I","A","N")

# Create dataframe of cell IDs and associated subtypes
subtypes_df <- list()
for(curr_subtype in names(clusters)){
  subtypes_df <- append(subtypes_df, list(cbind(colnames(all_samples_heatmap[,clusters[[curr_subtype]]]),curr_subtype)))
}
subtypes_df <- do.call(rbind,subtypes_df)
colnames(subtypes_df) <- c("cell_id","subtype")

# Add original order of cell in colData, then merge based on cell ID. Finally rearragne to orignal order
colData(ctcs)$original_order  <- 1:nrow(colData(ctcs))
subtype_order <- merge(colData(ctcs), subtypes_df, by="cell_id")
subtype_order <- subtype_order[order(subtype_order$original_order), ]

# Check that cell order is back to orginal order
all(subtype_order$cell_id == colData(ctcs)$cell_id)

# Add subtypes to colData
colData(ctcs)$subtype <- subtype_order$subtype

# Save data with subtype assignments
saveRDS(ctcs, "data/cytof_objects/ctcs_with_subtype.rds")
