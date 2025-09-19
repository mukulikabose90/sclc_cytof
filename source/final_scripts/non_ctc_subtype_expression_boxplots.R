################################################################################
# This script creates boxplots for expression of seleceted proteins in the 
# 4 "subtypes" detected in non CTCs
################################################################################
source("source/sclc_cytof_functions.R")

set.seed(42)
################################################################################
# Read in data
################################################################################
ctcs <- readRDS("data/cytof_objects/all_samples_non_ctcs_with_subtype.rds")

################################################################################
# Plot boxplots for SCLC TFs expression
################################################################################
markers_to_use <- c("NeuroD1","ASCL1","POU2F3","p-Rb")

p1 <- create_marker_boxplots(ctcs,markers_to_use,"subtype")

p1 <- p1+
  rremove("legend")+
  labs(x="Subtype",
       y="Expression")+
  scale_fill_manual(
    values = c("A" = "#dd4b33","N" = "#F1FAEE","P"= "#A8DADC", "I" = "#457B9D"))

################################################################################
# Plot boxplots for all proteins expression
################################################################################
markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam", "E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist", "SLUG", "PD-L1", "p-YAP", "CD44", "CD24")

p2 <- create_marker_boxplots(ctcs,markers_to_use,"subtype")

p2 <- p2+
  rremove("legend")+
  labs(x="Subtype",
       y="Expression")+
  scale_fill_manual(
    values = c("A" = "#dd4b33","N" = "#F1FAEE","P"= "#A8DADC", "I" = "#457B9D"))

################################################################################
# Save figures
################################################################################
tiff(glue("figures/non_ctcs_subtype_tf_expression_boxplot.tiff"), width=160,height=160, units = "mm", res=1000)
print(p1)
dev.off()

tiff(glue("figures/non_ctcs_subtype_expression_boxplot.tiff"), width=260,height=260, units = "mm", res=1000)
print(p2)
dev.off()

