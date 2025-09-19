# This script plots a heatmap of the SCLC subtype TFs expression. Then assigns each cell
# to a heatmap based on the expression of the TFs.

################################################################################
source("source/cytof_de_function.R")

col_fun = colorRamp2(c(-2, 0, 2), c("royalblue", "white", "firebrick"))
col_fun = colorRamp2(c(-2, 0, 2), c("royalblue4", "white", "firebrick4"))
col_fun = colorRamp2(c(-2, 0, 2), c("royalblue4","lightgray", "firebrick4"))

set.seed(42)
################################################################################
# Read in CyTOF data with cluster assignments
sce <- readRDS("data/cytof_objects/sclc_first_draw_with_clusters.rds")

# Subset to cancer cells in CTC cluster
# ctcs <- sce[,colData(sce)$new_clusters %in% c(1,6,10)]
ctcs <- sce[,colData(sce)$new_clusters %in% c(3,4)]
ctcs <- ctcs[,colData(ctcs)$condition == "cancer"]

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

rownames(heatmap_tf_scaled)



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
                                 show_annotation_name = T,
                                 annotation_legend_param = list(ncol=2))

# Create patient ID annotation
ttxx <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(heatmap_metadata$patient_id)))
names(ttxx) <- unique(heatmap_metadata$patient_id)
patient_anno <- HeatmapAnnotation("Patient ID" = heatmap_metadata$patient_id, 
                                  col = list("Patient ID"= ttxx),
                                  show_annotation_name = T,
                                  annotation_legend_param = list(ncol=3))


final_heatmap <- t(heatmap_tf_scaled)





naive_heatmap <- final_heatmap[,colnames(final_heatmap) %in% ctcs[,ctcs$treatment_status == "naive"]$cell_id]
treated_heatmap <- final_heatmap[,colnames(final_heatmap) %in% ctcs[,ctcs$treatment_status == "treated"]$cell_id]

# colnames(final_heatmap) <- NULL

ht <- Heatmap(final_heatmap, column_km = 4,top_annotation = patient_anno, name="Expression",
              cluster_columns = T, cluster_rows = F, show_column_names=F,col = col_fun,
              row_dend_reorder = F, column_title = c("SCLC-A", "SCLC-P","SCLC-I","SCLC-N"))

ht <- draw(ht)

ht <- Heatmap(naive_heatmap, column_km = 4, name="Expression",
              cluster_columns = T, cluster_rows = F, show_column_names=F,col = col_fun,
              row_dend_reorder = F)
naive_ht <- draw(ht)

ht <- Heatmap(treated_heatmap, column_km = 4, name="Expression",
              cluster_columns = T, cluster_rows = F, show_column_names=F,col = col_fun,
              row_dend_reorder = F)

treated_ht <- draw(ht)

length(column_order(naive_ht)$`1`)/ncol(naive_heatmap)
length(column_order(treated_ht)$`1`)/ncol(treated_heatmap)

a <- length(column_order(naive_ht)$`1`)
b <- ncol(naive_heatmap) - length(column_order(naive_ht)$`1`)
c <- length(column_order(treated_ht)$`1`)
d <- ncol(treated_heatmap) - length(column_order(treated_ht)$`1`)


contin <- matrix(c(a,c,b,d),ncol=2)

fisher.test(contin)

p1 <- fviz_nbclust(t(final_heatmap), kmeans, method='silhouette')+
  ggtitle("Optimal Number of Subclusters (Cancer Samples)")

png("figures/first_draw_ctcs_optimal_clusters.png", width = 12, height = 5, units = "in", res = 1200)
print(p1)
dev.off()

# Save figures
jpeg("figures/first_draw_ctcs_subtype_heatmap.jpg", width=200,height=100, units = "mm", res=1000)
print(ht)
dev.off()

# NAIVE SAMPLES
jpeg("figures/first_draw_ctcs_naive_subtype_heatmap.jpg", width = 12, height = 5, units = "in", res = 1200)
print(naive_ht)
dev.off()

#TREATED SAMPLES
jpeg("figures/first_draw_ctcs_treated_subtype_heatmap.jpg", width = 12, height = 5, units = "in", res = 1200)
print(treated_ht)
dev.off()

################################################################################
# Assign subtype to each cell based on kmeans clustering

# Get column row from heatmap
clusters <- column_order(ht)

# Rename clusters
names(clusters) <- c("A","P","I","N")

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
saveRDS(ctcs, "data/cytof_objects/first_draw_ctcs_with_subtype.rds")
