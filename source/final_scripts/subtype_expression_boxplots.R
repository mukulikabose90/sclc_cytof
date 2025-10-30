################################################################################
# This script creates boxplots for expression of seleceted proteins in the 
# 4 subtypes
################################################################################
source("source/sclc_cytof_functions.R")

set.seed(42)
################################################################################
# Read in data
################################################################################
ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

################################################################################
# Create boxplots for SCLC TFs in each subtype
################################################################################
markers_to_use <- c("NeuroD1","ASCL1","POU2F3","p-Rb")

p1 <- create_marker_boxplots(ctcs,markers_to_use,"subtype")

p1 <- p1+
  rremove("legend")+
  labs(x="Subtype",
       y="Expression")+
  scale_fill_manual(
    values = c("A" = "#dd4b33","N" = "#F1FAEE","P"= "#A8DADC", "Mes" = "#457B9D"))

################################################################################
# Create boxplots for selected proteins in each subtype
################################################################################
markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam", "E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist", "SLUG", "PD-L1", "p-YAP", "CD44", "CD24")

p2 <- create_marker_boxplots(ctcs,markers_to_use,"subtype")

p2 <- p2+
  rremove("legend")+
  labs(x="Subtype",
       y="Expression")+
  scale_fill_manual(
    values = c("A" = "#dd4b33","N" = "#F1FAEE","P"= "#A8DADC", "Mes" = "#457B9D"))


################################################################################
# Save figures
################################################################################
tiff(glue("figures/ctcs_subtype_tf_expression_boxplot.tiff"), width=160,height=160, units = "mm", res=1000)
print(p1)
dev.off()

tiff(glue("figures/ctcs_subtype_expression_boxplot.tiff"), width=260,height=260, units = "mm", res=1000)
print(p2)
dev.off()

 