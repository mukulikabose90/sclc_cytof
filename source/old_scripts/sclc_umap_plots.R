# This script plots the UMAP plots for the SCLC CyTOF data.
sce <- readRDS("data/cytof_objects/sclc_cytof_sce_object.rds")

sce <- sce[,colData(sce)$sample_type != "cell_line"]




# sce <- sce[,colData(sce)$condition == "cancer"]
# PLOT UMAP PLOTS

set.seed(42)
sce <- runDR(sce, "UMAP", cells = 5e3, features = "state")

plotDR(sce, "UMAP", color_by = "patient_id", scale = T)+
  geom_point(size=1)

plotDR(sce, "UMAP", color_by = "condition", scale = T)+
  geom_point(size=1)

plotDR(sce, "UMAP", color_by = "NeuroD1",facet_by ="condition", scale = T)+
  geom_point(size=1)


plotDR(sce, "UMAP", color_by = "patient_id", facet_by="patient_id",scale = T)+
  geom_point(size=1)
