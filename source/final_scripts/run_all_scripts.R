renv::restore()

if (!dir.exists("figures/")) {
  dir.create("figures/")
}

if (!dir.exists("data/")) {
  dir.create("data/")
}
################################################################################
# Process and QC
################################################################################
source("source/final_scripts/process_samples.R")
source("source/final_scripts/cytof_samples_ qc.R")

################################################################################
# Identify cancer-enriched populations
################################################################################
source("source/final_scripts/all_samples_clustering.R")
source("source/final_scripts/differential_abundance_lm.R")
source("source/final_scripts/all_samples_clustering_opacity.R")

source("source/final_scripts/all_samples_cluster_expression_heatmap.R")
source("source/final_scripts/all_samples_cluster_expression_boxplots.R")

################################################################################
# CTC detection
################################################################################
source("source/final_scripts/ce_subclustering.R")
source("source/final_scripts/ce_subcluster_expression_heatmap.R")
source("source/final_scripts/ce_subcluster_expression_boxplots.R")

################################################################################
# pRb of cell groups
################################################################################
source("source/final_scripts/ctcs_vs_normal_cells_expression.R")

################################################################################
# CTC analysis
################################################################################
source("source/final_scripts/ctcs_subtype_detection_heatmap.R")
source("source/final_scripts/subtype_expression_heatmap.R")
source("source/final_scripts/subtype_expression_boxplots.R")

source("source/final_scripts/non_ctc_subtype_detection_heatmap.R")
source("source/final_scripts/non_ctc_subtype_expression_heatmap.R")
source("source/final_scripts/non_ctc_subtype_expression_boxplots.R")

source("source/final_scripts/normal_subtype_detection_heatmap.R")
source("source/final_scripts/normal_subtype_expression_heatmap.R")
source("source/final_scripts/normal_subtype_expression_boxplots.R")

source("source/final_scripts/all_cells_subtype_treatment_status_three_comparison.R")

################################################################################
# Treatment status analysis
################################################################################
source("source/final_scripts/treatment_status_expr_boxplots.R")
source("source/final_scripts/treatment_status_expr_heatmap.R")
source("source/final_scripts/treatment_status_subtype_barplots.R")

################################################################################
# Tarla analysis
################################################################################
source("source/final_scripts/pre_post_tarla_expr_boxplots.R")
source("source/final_scripts/pre_post_tarla_expr_heatmap.R")
source("source/final_scripts/pre_post_tarla_subtype_barplots.R")

################################################################################
# Longitudinal analysis
################################################################################
source("source/final_scripts/longitudinal_subtype_barplots.R")
source("source/final_scripts/longitudinal_expression_heatmap.R")

source("source/final_scripts/ctc_subtype_treatment_proportion_barplot.R")

 ################################################################################
# Alluvial Plots
################################################################################
source("source/final_scripts/subtype_treatment_alluvial.R")
source("source/final_scripts/downsample_alluvial.R")

################################################################################
# Downsampling analysis
################################################################################
source("source/final_scripts/bootstrap_subtype_sampling.R")
source("source/final_scripts/downsampled_subtype_treatment_status_three_comparison.R")
source("source/final_scripts/downsampled_samples_subtype_barplots.R")

################################################################################
# Xenograft samples plots
################################################################################
source("source/final_scripts/xenograft_samples_barplots.R")
source("source/final_scripts/xenograft_samples_violin_plots.R")
source("source/final_scripts/xenograft_samples_heatmap.R")

source("source/final_scripts/pleural_effusion_heatmap.R")

################################################################################
# Calculate percent CTCs
################################################################################
source("source/final_scripts/calc_percent_ctcs.R")



