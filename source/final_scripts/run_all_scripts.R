# Process and QC
source("source/final_scripts/process_samples.R") 
source("source/final_scripts/cytof_samples_qc.R")

# Identify cancer-enriched populations
source("source/final_scripts/all_samples_clustering.R")
source("source/final_scripts/cluster_condition_chi_square_plot.R")
source("source/final_scripts/all_samples_clustering_opacity.R")

source("source/final_scripts/all_samples_cluster_expression_heatmap.R")
source("source/final_scripts/all_samples_cluster_expression_boxplots.R")

# CTC detection 
source("source/final_scripts/ce_subclustering.R") 
source("source/final_scripts/ce_subcluster_expression_heatmap.R")
source("source/final_scripts/ce_subcluster_expression_boxplots.R")
source("source/final_scripts/ce_expression_umaps.R")

# CTC analysis
source("source/final_scripts/ctcs_subtype_detection_heatmap.R")
source("source/final_scripts/subtype_expression_heatmap.R")
source("source/final_scripts/subtype_expression_boxplots.R")

source("source/final_scripts/non_ctc_subtype_detection_heatmap.R")
source("source/final_scripts/non_ctc_subtype_expression_heatmap.R")
source("source/final_scripts/non_ctc_subtype_expression_boxplots.R")

source("source/final_scripts/normal_subtype_detection_heatmap.R")
source("source/final_scripts/normal_subtype_expression_heatmap.R")
source("source/final_scripts/normal_subtype_expression_boxplots.R")

source("source/final_scripts/subtype_treatment_association.R")

# Treatment status analysis 
source("source/final_scripts/treatment_status_expr_boxplots.R")
source("source/final_scripts/treatment_status_expr_heatmap.R")
source("source/final_scripts/treatment_status_subtype_barplots.R")

#Tarla analysis
source("source/final_scripts/pre_post_tarla_expr_boxplots.R")
source("source/final_scripts/pre_post_tarla_expr_heatmap.R")
source("source/final_scripts/pre_post_tarla_subtype_barplots.R")


#Longitudinal analysis
source("source/final_scripts/longitudinal_subtype_barplots.R")
source("source/final_scripts/longitudinal_expression_heatmap.R")


# Protein Expression UMAPs
# source("source/final_scripts/all_samples_ctc_expression_umaps.R")
# source("source/final_scripts/tarla_ctc_expression_umaps.R")



