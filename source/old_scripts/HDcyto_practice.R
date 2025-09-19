library(readxl)
library(HDCytoData)
library(CATALYST)
library(diffcyt)

################################################################################

# Download Bodenmiller metadata
url <- "https://zenodo.org/records/10039274/files"
md <- "PBMC8_metadata.xlsx"
download.file(file.path(url, md), destfile = md, mode = "wb")
md <- read_excel(md)
head(data.frame(md))

# Load Bodenmiller flow dataset
fs <- Bodenmiller_BCR_XL_flowSet()

# Read in panel data
panel <- "PBMC8_panel_v3.xlsx"
download.file(file.path(url, panel), destfile = panel, mode = "wb")
panel <- read_excel(panel)
head(data.frame(panel))


fs[[1]]
pData(parameters(fs[[1]]))
colnames(fs)
all(panel$fcs_colname %in% colnames(fs))


# specify levels for conditions & sample IDs to assure desired ordering
md$condition <- factor(md$condition, levels = c("Ref", "BCRXL"))
md$sample_id <- factor(md$sample_id, 
                       levels = md$sample_id[order(md$condition)])

# construct SingleCellExperiment
sce <- prepData(fs, panel, md, features = panel$fcs_colname)

sce@metadata$experiment_info

rownames(sce@assays@data$counts)

p <- plotExprs(sce, color_by = "condition")
p$facet$params$ncol <- 6
# p

n_cells(sce) 

# plotCounts(sce, group_by = "sample_id", color_by = "condition")

# pbMDS(sce, color_by = "condition", label_by = "sample_id")

# plotExprHeatmap(sce, scale = "last",
#                 hm_pal = rev(hcl.colors(10, "YlGnBu")))


# plotNRS(sce, features = "type", color_by = "condition")


set.seed(1234)
sce <- cluster(sce, features = "type",
               xdim = 10, ydim = 10, maxK = 20, seed = 1234)


# plotExprHeatmap(sce, features = "type", 
#                 by = "cluster_id", k = "meta20", 
#                 bars = TRUE, perc = TRUE)

# plotClusterExprs(sce, k = "meta20", features = "type")



# run t-SNE/UMAP on at most 500/1000 cells per sample
set.seed(1234)
sce <- runDR(sce, "UMAP", cells = 1e3, features = "type")

plotDR(sce, "UMAP", color_by = "CD4")


# Merging clusters

merging_table1 <- "PBMC8_cluster_merging1.xlsx"
download.file(file.path(url, merging_table1), 
              destfile = merging_table1, mode = "wb")
merging_table1 <- read_excel(merging_table1)
head(data.frame(merging_table1))
# convert to factor with merged clusters in desired order
merging_table1$new_cluster <- factor(merging_table1$new_cluster, 
                                     levels = c("B-cells IgM+", "B-cells IgM-", "CD4 T-cells",
                                                "CD8 T-cells", "DC", "NK cells", "monocytes", "surface-"))

# apply manual merging
sce <- mergeClusters(sce, k = "meta20", 
                     table = merging_table1, id = "merging1")

plotDR(sce, "UMAP", color_by = "merging1")

plotAbundances(sce, k = "merging1", by = "sample_id")

plotAbundances(sce, k = "merging1", by = "cluster_id", shape_by = "patient_id")


ei <- metadata(sce)$experiment_info
(da_formula1 <- createFormula(ei, 
                              cols_fixed = "condition", 
                              cols_random = "sample_id"))


contrast <- createContrast(c(0, 1))


da_res1 <- diffcyt(sce, 
                   formula = da_formula1, contrast = contrast,
                   analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",
                   clustering_to_use = "merging1", verbose = FALSE)


rowData(da_res1$res) 

og_ei <- ei
og_form <- da_formula1
og_contrast <- contrast
og_res <- da_res1

FDR_cutoff <- 0.05

plotDiffHeatmap(sce, rowData(da_res1$res), all = TRUE, fdr = FDR_cutoff)



