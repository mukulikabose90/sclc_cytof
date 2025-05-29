
sce <- readRDS("data/cytof_objects/sclc_all_samples_object_no_qc.rds")

################################################################################
# Remove cell line samples
################################################################################
blood_samples <- as.data.frame(sce@colData) %>%
  dplyr::filter(sample_type == "blood") %>%
  pull(collection_id) %>%
  as.character()

blood_samples <- as.data.frame(sce@colData) %>%
  dplyr::filter(sample_type == "cell_line") %>%
  pull(collection_id) %>%
  as.character()

sce <- sce[,sce$collection_id %in% blood_samples]

################################################################################

################################################################################
# markers <- as.data.frame(rowData(sce)) %>%
#   dplyr::filter(marker_class == "state") %>%
#   pull(marker_name)
# 
# temp <- sce[markers,sce$collection_id == "H1105-1"]
# 
# p1 <- plotExprs(temp, color_by = "experiment_id", assay = "exprs")
# 
# p1

# temp <- sce[markers,sce$collection_id == "NJH29-1"]
# 
# p1 <- plotExprs(temp, color_by = "experiment_id", assay = "exprs")
# # 
# temp <- sce_corrected[markers,sce_corrected$collection_id == "NJH29-1"]
# 
# p2 <- plotExprs(temp, color_by = "experiment_id", assay = "exprs")
# 
# p1+ggtitle("NJH29 (no batch correction)")
# p2+ggtitle("NJH29 (batch corrected)")
# 
# ################################################################################

# Get samples that are run in multiple experiments
# to_test <- as.data.frame(sce@colData) %>%
#   dplyr::select(experiment_id,collection_id) %>%
#   distinct() %>%
# dplyr::count(collection_id) %>%
#   arrange(desc(n)) %>%
#   dplyr::filter(n > 1) %>%
#   pull(collection_id) %>%
#   as.character()
# 
# 
# sort(to_test)
# 
# # Checking SC454-1
# temp <- sce[markers,sce$collection_id == "SC454-1"]
# 
# plotExprs(temp, color_by = "experiment_id", assay = "exprs")
# 
# temp <- sce_corrected[markers,sce_corrected$collection_id == "SC414-1"]
# plotExprs(temp, color_by = "experiment_id",assay = "exprs")
################################################################################
# Proportional downsampling
################################################################################
downsample_df <- as.data.frame(colData(sce)) %>% 
  dplyr::count(collection_id,experiment_id) %>% 
  dplyr::count(collection_id) %>% 
  arrange(desc(n))

downsampled_objs <- list()
for(i in 1:nrow(downsample_df)){
  
  curr_sce <- sce[,sce$collection_id == downsample_df[i,1]]
  
  downsample_frac <- 1/downsample_df[i,2]
  
  max_n <- downsample_frac*ncol(curr_sce)
  
  curr_sce <- downsampleSCE(curr_sce,
                            group_by = "collection_id",
                            maxN = max_n)
  
  downsampled_objs <- append(downsampled_objs,list(curr_sce))
}

sce <- do.call(cbind, downsampled_objs)

################################################################################
# markers <- as.data.frame(rowData(sce)) %>%
#   dplyr::filter(marker_class == "state") %>%
#   pull(marker_name)
# 
# temp <- sce[markers,sce$collection_id == "H1105-1"]
# 
# p1 <- plotExprs(temp, color_by = "experiment_id", assay = "exprs")
# 
# p1

# temp <- sce[markers,sce$collection_id == "NJH29-1"]
# 
# p1 <- plotExprs(temp, color_by = "experiment_id", assay = "exprs")
# # 
# temp <- sce_corrected[markers,sce_corrected$collection_id == "NJH29-1"]
# 
# p2 <- plotExprs(temp, color_by = "experiment_id", assay = "exprs")
# 
# p1+ggtitle("NJH29 (no batch correction)")
# p2+ggtitle("NJH29 (batch corrected)")
# 
# ################################################################################

# Get samples that are run in multiple experiments
# to_test <- as.data.frame(sce@colData) %>%
#   dplyr::select(experiment_id,collection_id) %>%
#   distinct() %>%
# dplyr::count(collection_id) %>%
#   arrange(desc(n)) %>%
#   dplyr::filter(n > 1) %>%
#   pull(collection_id) %>%
#   as.character()
# 
# 
# sort(to_test)
# 
# # Checking SC454-1
# temp <- sce[markers,sce$collection_id == "SC454-1"]
# 
# plotExprs(temp, color_by = "experiment_id", assay = "exprs")
# 
# temp <- sce_corrected[markers,sce_corrected$collection_id == "SC414-1"]
# plotExprs(temp, color_by = "experiment_id",assay = "exprs")

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

covariates <- c("condition","experiment_id")
covariates <- c("experiment_id")
covariates <- paste(covariates, collapse = " + ")

all_markers <- markers_to_use

pre_models <- list()
for(curr_marker in all_markers){
  
  expr_mat <- sce[curr_marker,]@assays@data$exprs
  
  rownames(expr_mat) <- "curr_marker"
  df <- cbind(as.data.frame(colData(sce)),t(expr_mat))
  
  res <- lm(formula = as.formula(paste0("curr_marker ~ ", covariates)), data = df)
  
  pre_models <- append(pre_models, list((res)))
}


for(i in 1:length(pre_models)){
  print(pre_models[[i]])
  cat("\n \n")
}

################################################################################
# cyCombine batch correction
################################################################################
marker_info <- read.csv("data/cytof_panel_info.csv")
marker_info <- data.frame(marker_info, stringsAsFactors = FALSE)

markers <- marker_info %>%
  dplyr::filter(marker_class != "none") %>%
  pull(antigen)


y <- assay(sce, "exprs")

colnames(y) <- paste0("cell_",1:ncol(y))

df <- data.frame(t(y), colData(sce), check.names = FALSE)

colnames(df)[which(colnames(df) == "experiment_id")] <- "batch"

corrected <- batch_correct(df,
                           markers = markers)

corrected <- as.matrix(t(corrected[,3:40]))

corrected <- corrected[rownames(y),]

assay(sce, "exprs") <- corrected

################################################################################

covariates <- c("condition","experiment_id")
covariates <- c("experiment_id")
covariates <- paste(covariates, collapse = " + ")

all_markers <- markers_to_use

post_models <- list()
for(curr_marker in all_markers){
  
  expr_mat <- sce[curr_marker,]@assays@data$exprs
  
  rownames(expr_mat) <- "curr_marker"
  df <- cbind(as.data.frame(colData(sce)),t(expr_mat))
  
  res <- lm(formula = as.formula(paste0("curr_marker ~ ", covariates)), data = df)
  
  post_models <- append(post_models, list((res)))
}





get_mean_abs_batch_coef <- function(model) {
  coefs <- coef(model)
  batch_coefs <- coefs[grepl("^experiment_id", names(coefs))]
  return(mean(abs(batch_coefs)))
}

final_df <- list()
for(i in 1:length(all_markers)){
  
  curr_marker <- all_markers[i]
  
  pre_mean <- get_mean_abs_batch_coef(pre_models[[i]])
  post_mean <- get_mean_abs_batch_coef(post_models[[i]])
  
  
  final_df[["marker"]] <- append(final_df[["marker"]], rep(curr_marker,2))
  final_df[["class"]] <- append(final_df[["class"]], c("pre","post"))
  final_df[["value"]] <- append(final_df[["value"]], c(pre_mean,post_mean))
  
}

final_df <- as.data.frame(final_df)

final_df$class <- factor(final_df$class, levels = c("pre","post"))

ggplot(final_df)+
  geom_col(aes(x=marker,y=value,fill=class), position = "dodge")+
  labs(x="Protein",
       y="Mean |Coefficient|",
       fill = "")









