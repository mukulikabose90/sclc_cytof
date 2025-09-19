# This script calculates the median expression for each SCLC subtype TF for each
# patient. It then plots a scaled heatmap.

################################################################################
source("source/cytof_de_function.R")

col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))

set.seed(42)
################################################################################
# Read in CyTOF data with cluster assignments
sce <- readRDS("data/cytof_objects/sclc_cytof_sce_object_with_clusters.rds")

# Subset to cell in cluster 1
ctcs <- sce[,colData(sce)$new_clusters == 1]

ctcs <- ctcs[,colData(ctcs)$condition == "cancer"]

# Subset to only state markers
ctcs <- ctcs[rowData(ctcs)$marker_class == "state",]


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



# Create heatmap of median expression for each patient
heatmap <- gg_df %>% 
  group_by(patient_id, antigen) %>% 
  summarise(mean(expression)) %>% 
  pivot_wider(names_from = "antigen", values_from = "mean(expression)") %>% 
  column_to_rownames("patient_id")

# heatmap_scaled <-  apply(heatmap, MARGIN = 2, FUN = function(x) scale_values(x))
heatmap_scaled <- scale(heatmap)

# Subset to only SCLC subtype TFs and cancer samples

# heatmap_tf <- heatmap[,colnames(heatmap) %in% c("POU2F3","NeuroD1","p-YAP","ASCL1")]
# heatmap_tf <- heatmap[,colnames(heatmap) %in% c("POU2F3","NeuroD1","ASCL1")]
# heatmap_tf <- heatmap_tf[grepl("SC",rownames(heatmap_tf)),]

# Scale expression for each protein between 0 and 1
# heatmap_tf_scaled <- apply(heatmap_tf, MARGIN = 2, FUN = function(x) scale_values(x))

heatmap_tf_scaled <- heatmap_scaled[,colnames(heatmap_scaled) %in% c("POU2F3","NeuroD1","ASCL1")]
##################################################################



final_heatmap <- t(heatmap_tf_scaled)
Heatmap(final_heatmap, column_km = 4, name="Scaled\nExpression",
        cluster_columns = T, cluster_rows = F)





heatmap_tf_scaled[1,]

hist(as.vector(heatmap_tf_scaled))


binary_matrix <- ifelse(heatmap_tf_scaled > .5, 1, 0)
Heatmap(binary_matrix, name = "Median\nExpression", col = col_fun,
        cluster_columns = F, cluster_rows = T)





subtypes <- c()
for(i in 1:nrow(binary_matrix)){
  value <- colnames(heatmap_tf)[which(binary_matrix[i,] == 1)]
  
  subtypes[i] <- ifelse(identical(value, character(0)), "I",value)
  
}


temp <- merge(colData(sce),cbind("patient_id"=rownames(binary_matrix),subtypes), by="patient_id", all=T)

colData(sce)$subtypes <- temp$subtypes

saveRDS(sce, "data/cytof_objects/sclc_cytof_sce_object_with_subtypes.rds")


