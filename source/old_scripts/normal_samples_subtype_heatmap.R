# This script plots a heatmap of the SCLC subtype TFs expression. Then assigns each cell
# to a heatmap based on the expression of the TFs.

################################################################################
source("source/cytof_de_function.R")

col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))

set.seed(42)
################################################################################
# Read in CyTOF data with cluster assignments
sce <- readRDS("data/cytof_objects/sclc_cytof_sce_object_with_clusters.rds")

# Subset to cancer cells in CTC cluster
ctcs <- sce[,!colData(sce)$new_clusters %in% c(4,5,6,7)]
# ctcs <- sce[,colData(sce)$new_clusters %in% c(2,3,6,7)]
ctcs <- ctcs[,colData(ctcs)$condition == "normal"]



# Subset to only state markers
ctcs <- ctcs[rowData(ctcs)$marker_class == "state",]

# Give each cell an ID
colData(ctcs)$cell_id <- 1:nrow(colData(ctcs))
colData(ctcs)$cell_id <- paste0("cell_",1:nrow(colData(ctcs)))

# Remove patients with less than 10 cells
patients_to_keep <- as.data.frame(colData(ctcs)) %>% 
  dplyr::count(patient_id) %>% 
  dplyr::filter(n > 10) %>%
  pull(patient_id) %>% 
  as.character()

ctcs <- ctcs[,colData(ctcs)$patient_id %in% patients_to_keep]

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
  dplyr::select(cell_id,collection_id,sample_id,patient_id,condition,sample_num, treatment_status) %>% 
  select(-c(sample_id)) %>% 
  distinct()

# Create collection ID annotation
ttxx <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(heatmap_metadata$collection_id)))
names(ttxx) <- unique(heatmap_metadata$collection_id)
sample_anno <- HeatmapAnnotation("Sample ID" = heatmap_metadata$collection_id, 
                                 col = list("Sample ID"= ttxx),
                                 show_annotation_name = T)

# Create patient ID annotation
ttxx <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(heatmap_metadata$patient_id)))
names(ttxx) <- unique(heatmap_metadata$patient_id)
patient_anno <- HeatmapAnnotation("Patient ID" = heatmap_metadata$patient_id, 
                                  col = list("Patient ID"= ttxx),
                                  show_annotation_name = T)


final_heatmap <- t(heatmap_tf_scaled)
# colnames(final_heatmap) <- NULL



ht <- Heatmap(final_heatmap,column_km=4,top_annotation = patient_anno, name="Expression",
              cluster_columns = T, cluster_rows = F, show_column_names=F,
              row_dend_reorder = F)

ht <- draw(ht)

p1 <- fviz_nbclust(t(final_heatmap), kmeans, method='silhouette')+
  ggtitle("Optimal Number of Subclusters (Normal Samples)")

png("figures/normal_optimal_clusters.png", width = 12, height = 5, units = "in", res = 1200)
print(p1)
dev.off()

# Save figures
png("figures/normals_subtype_heatmap.png", width = 12, height = 5, units = "in", res = 1200)
# tiff("figures/ctcs_subtype_heatmap.tiff", width = 12, height = 5, units = "in", res = 1200)
draw(ht)
dev.off()
# dev.off()

################################################################################
# Assign subtype to each cell based on kmeans clustering

# Get column row from heatmap
clusters <- column_order(ht)

names(column_order(ht)) <- c("SCLC-P","SCLC-I","SCLC-A","SCLC-N")

# Rename clusters
names(clusters) <- c("P","I","A","N")

# Create dataframe of cell IDs and associated subtypes
subtypes_df <- list()
for(curr_subtype in names(clusters)){
  subtypes_df <- append(subtypes_df, list(cbind(colnames(final_heatmap[,clusters[[curr_subtype]]]),curr_subtype)))
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
saveRDS(ctcs, "data/cytof_objects/normals_with_subtype.rds")
