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

colData(ctcs)$cell_id <- 1:nrow(colData(ctcs))

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
  select(cell_id,antigen,expression) %>% 
  pivot_wider(names_from = antigen, values_from = expression) %>% 
  column_to_rownames("cell_id")

# rownames(heatmap) <- NULL

dim(heatmap)

# heatmap_scaled <-  apply(heatmap, MARGIN = 2, FUN = function(x) scale_values(x))

heatmap_scaled <- scale(heatmap)

# Subset to only SCLC subtype TFs and cancer samples


# Scale expression for each protein between 0 and 1
# heatmap_tf_scaled <- apply(heatmap_tf, MARGIN = 2, FUN = function(x) scale_values(x))
# heatmap_tf <- heatmap[,colnames(heatmap) %in% c("POU2F3","NeuroD1","ASCL1")]
# heatmap_tf_scaled <- scale(heatmap_tf)
heatmap_tf_scaled <- heatmap_scaled[,colnames(heatmap_scaled) %in% c("POU2F3","NeuroD1","ASCL1")]


heatmap_metadata <- data.frame(colData(ctcs)) %>% 
  dplyr::select(cell_id,patient_id) %>% 
  distinct()

ttxx <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(heatmap_metadata$patient_id)))
names(ttxx) <- unique(heatmap_metadata$patient_id)
patient_anno <- HeatmapAnnotation("Patient ID" = heatmap_metadata$patient_id, 
                              col = list("Patient ID"= ttxx),
                              show_annotation_name = T)


final_heatmap <- t(heatmap_tf_scaled)
colnames(final_heatmap) <- NULL
ht <- Heatmap(final_heatmap, column_km = 4,top_annotation = patient_anno, name="Scaled\nExpression",
        cluster_columns = T, cluster_rows = F)


ht <- draw(ht)

clusters <- column_order(ht)

names(clusters) <- c("P","I","A","N")

subtypes_df <- list()
for(curr_subtype in names(clusters)){
  
  subtypes_df <- append(subtypes_df, list(cbind(clusters[[curr_subtype]],curr_subtype)))
  
}

subtypes_df <- do.call(rbind,subtypes_df)
colnames(subtypes_df) <- c("cell_id","subtype")



subtype_order <- merge(colData(ctcs), subtypes_df, by="cell_id",all.x=T)

colData(ctcs)$subtype <- subtype_order$subtype

temp <- as.data.frame(colData(ctcs)) %>% 
  dplyr::count(patient_id,subtype) %>% 
  dplyr::filter(!is.na(subtype))


ggplot(temp)+
  geom_col(aes(x=patient_id,y=n,fill=subtype), position = "dodge")


cluster_prop_df <- as.data.frame(colData(ctcs)) %>%
  group_by(patient_id,subtype) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))%>% 
  group_by(patient_id) %>% 
  mutate(total = sum(n))


ggplot(cluster_prop_df,aes(x=patient_id,y=freq, fill=subtype))+
  geom_col()+
  geom_text(aes(y = 1.03,label=total))+
  xlab("Patient")+
  ylab("Percentage")+
  theme(axis.text.x = element_text(angle=70,hjust=1,vjust=1))


cluster_prop_df <- as.data.frame(colData(ctcs)) %>%
  group_by(stage,subtype) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))%>% 
  group_by(stage) %>% 
  mutate(total = sum(n))




ggplot(cluster_prop_df,aes(x=stage,y=freq, fill=subtype))+
  geom_col(position = "dodge")+
  geom_text(aes(y = 1.03,label=total))+
  xlab("Patient")+
  ylab("Percentage")+
  theme(axis.text.x = element_text(angle=70,hjust=1,vjust=1))


ggplot(cluster_prop_df) +
  geom_boxplot(aes(x=subtype, y=freq))
  facet_wrap(~subtype)

