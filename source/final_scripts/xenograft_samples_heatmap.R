################################################################################
# This script plots a heatmap for protein expression for the patients that had
# PDXs 
################################################################################
source("source/sclc_cytof_functions.R")

set.seed(42)
################################################################################
# Read in data
################################################################################
ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

################################################################################
# Create heatmap
################################################################################
markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam", "E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist", "SLUG", "PD-L1", "p-YAP", "CD44", "CD24")

col_fun = colorRamp2(c(-2, -1, 0, 1, 2), 
                     c("#313695",  # deep blue
                       "#74add1",  # light blue
                       "#f7f7f7",  # white (center, 0)
                       "#f46d43",  # light red
                       "#a50026"))


patients_to_use <- c("MDA-SC293","MDA-SC506","MDA-SC443")

curr_data <- ctcs

#Scale expression
assay(curr_data, "exprs") <- t(scale(t(assay(curr_data, "exprs"))))

curr_data <- curr_data[,curr_data$patient_id %in% patients_to_use & curr_data$treatment_status == "naive"]

curr_heatmap <- matrix(NA,ncol=3,nrow=length(markers_to_use))
for(curr_patient in patients_to_use){
  
  i <- which(curr_patient == patients_to_use)
    
  y <- assay(curr_data[,curr_data$patient_id == curr_patient], "exprs")
  
  y <- y[markers_to_use,]
  
  curr_heatmap[,i] <- rowMedians(y)
  
  
}

colnames(curr_heatmap) <- patients_to_use
rownames(curr_heatmap) <- markers_to_use

################################################################################
# Plot heatmap
################################################################################
curr_ht <- Heatmap(curr_heatmap, cluster_rows = F, cluster_columns = F,column_names_rot = 0,column_names_centered = T,
                   name="Median Scaled\n   Expression", column_title = "", col = col_fun)


p1 <- draw(curr_ht)

################################################################################
# Save figure
################################################################################
tiff(glue("figures/xenograft_expression_heatmap.tiff"), width=160,height=100, units = "mm", res=600)
print(p1)
dev.off()
