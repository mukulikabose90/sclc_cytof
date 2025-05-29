

sce <- readRDS("data/cytof_objects/sce_not_sampled.rds")

################################################################################
# Remove outlier experiments
################################################################################
sce <- sce[,sce$experiment_id != "531050"]
# sce <- sce[,sce$experiment_id != "508814"]
# sce <- sce[,sce$experiment_id != "513549"]

################################################################################
# Remove low viability samples
################################################################################
low_viability_samples <- paste0("NORMAL", c(10,11,15,18))

sce <- sce[,!sce$patient_id %in% low_viability_samples]

################################################################################
# Remove blood bank samples
################################################################################
# blood_bank_samples <- paste0("NORMAL", 7:20)
# 
# length(unique(sce$patient_id))
# 
# sce <- sce[,!sce$patient_id %in% blood_bank_samples]
# 
# as.character(unique(colData(sce)$patient_id))

################################################################################
# remove collections with < 30 cells
################################################################################
samples_to_remove <- as.data.frame(sce@colData) %>%
  dplyr::count(collection_id) %>%
  dplyr::filter(n<30) %>%
  pull(collection_id) %>%
  as.character()

sce <- sce[,!sce$collection_id %in% samples_to_remove]

################################################################################
# cyCombine batch correction
################################################################################
marker_info <- read.csv("data/cytof_panel_info.csv")
marker_info <- data.frame(marker_info, stringsAsFactors = FALSE)

markers <- marker_info %>%
  dplyr::filter(marker_class != "none") %>%
  pull(antigen)

# get expression data
y <- assay(sce, "exprs")

colnames(y) <- paste0("cell_",1:ncol(y))

df <- data.frame(t(y), colData(sce), check.names = FALSE)

colnames(df)[which(colnames(df) == "experiment_id")] <- "batch"

corrected <- batch_correct(df,markers = markers)

corrected <- as.matrix(t(corrected[,3:40]))

corrected <- corrected[rownames(y),]

assay(sce, "exprs") <- corrected



sce <- CATALYST::cluster(sce, features = "state",
                         xdim = 10, ydim = 10, maxK = 20, seed = script_seed)

sce <- runDR(sce, "UMAP", cells = 5e3, features = "state")

sce@metadata$delta_area





sce_query <- sce
sce_ref <- readRDS("data/cytof_objects/sclc_all_samples_with_clusters.rds")


# Extract expression matrix (transpose: cells x features)
ref_exprs <- t(assay(sce_ref, "exprs"))  # replace "exprs" with appropriate assay

# Run PCA using prcomp (centered and scaled)
pca_ref_model <- prcomp(ref_exprs, center = TRUE, scale. = TRUE)


# Project reference cells into PCA space
pca_ref <- pca_ref_model$x[, 1:10]  # top 10 PCs

query_exprs <- t(assay(sce_query, "exprs"))

# Center and scale query data using reference PCA parameters
query_scaled <- scale(query_exprs,
                      center = pca_ref_model$center,
                      scale = pca_ref_model$scale)

# Project query cells into the same PCA space
pca_query <- as.matrix(query_scaled) %*% pca_ref_model$rotation[, 1:10]

library(FNN)

# Find nearest neighbor in reference PCA space
nn <- get.knnx(pca_ref, pca_query, k = 1)

# Map to reference clusters
ref_clusters <- colData(sce_ref)$new_clusters  # replace with your actual column name
colData(sce_query)$projected_cluster <- ref_clusters[nn$nn.index[, 1]]


# Store PCA
reducedDim(sce_query, "PCA") <- pca_query

# Run UMAP on projected PCA
library(scater)
sce_query <- runUMAP(sce_query, dimred = "PCA")

plotUMAP(sce_query, color_by="projected_cluster")
plotUMAP(sce_ref, color_by="new_clusters")



sce_query %>% 
  colData() %>% 
  as.data.frame() %>% 
  filter(projected_cluster %in% ctc_clusters) %>% 
  count(condition)



sce_query[,sce_query$projected_cluster %in% ctc_clusters & sce_query$condition == "cancer"]



sce <- sce_query

ctc_clusters <- readRDS("data/ctc_clusters.rds")

sce$ctc <- ifelse(sce$projected_cluster %in% ctc_clusters & sce$condition == "cancer", "ctc", "non-ctc")
# sce$ctc <- ifelse(sce$new_clusters %in% ctc_clusters, "ctc", "non-ctc")

y <- assay(sce, "exprs")

df <- data.frame(t(y), colData(sce), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")

gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(sce)))

sclc_tfs <- c("NeuroD1","ASCL1","POU2F3","p-Rb")

sclc_tfs <- markers

temp <- gg_df %>%
  dplyr::filter(antigen %in% sclc_tfs)

temp$ctc <- factor(temp$ctc, levels=c("ctc","non-ctc"))


p <- ggboxplot(temp, x="projected_cluster",y="expression", fill="ctc")

p1 <- p + facet_wrap(~antigen, scales = "free_y")+
  stat_compare_means()+
  scale_fill_manual(values=c("#E57373","#64B5F6"))

p <- ggboxplot(temp, x="ctc",y="expression", fill="ctc")

p2 <- p + facet_wrap(~antigen, scales = "free_y")+
  stat_compare_means()+
  scale_fill_manual(values=c("#E57373","#64B5F6"))


p1

p2

