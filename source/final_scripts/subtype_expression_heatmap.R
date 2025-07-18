source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

table(ctcs$subtype)

# assay(ctcs, "exprs") <- t(scale(t(assay(ctcs, "exprs"))))

markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam", "E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist", "SLUG", "PD-L1", "p-YAP", "CD44", "CD24")

ht <- create_expression_heatmap(ctcs, "subtype", markers_to_use,"", scale = T)

print(ht)

jpeg("figures/ctcs_subtype_expression_heatmap.jpg", width=200,height=100, units = "mm", res=1000)
print(ht)
dev.off()

# 
# # assay(ctcs, "exprs") <- t(scale(t(assay(ctcs, "exprs"))))
# 
# all_markers <- readRDS("data/state_markers.rds")
# 
# naive <- ctcs[,ctcs$treatment_status == "naive"]
# naive_ht <- create_expression_heatmap(naive, "subtype", all_markers)
# 
# treated <- ctcs[,ctcs$treatment_status == "treated"]
# treated_ht <- create_expression_heatmap(treated, "subtype", all_markers)
# 
# naive_ht+treated_ht
# 
# subtypes <- c("A","N","P","I")
# heatmaps <- list()
# for(curr_subtype in subtypes){
#   curr_sce <- ctcs[,ctcs$subtype == curr_subtype]
#   
#   # curr_sce <- curr_sce[-which(rownames(curr_sce) %in% c("CD24","Vimentin"))]
#   
#   curr_ht <- create_expression_heatmap(curr_sce, "treatment_status", rownames(curr_sce),curr_subtype,
#                                        scale = F)
#  
#   heatmaps <- append(heatmaps, list(curr_ht))
#    
# }
# 
# 
# 
# heatmaps[[1]]+heatmaps[[2]]+heatmaps[[3]]+heatmaps[[4]]
# # 
# # jpeg("figures/ctcs_subtype_expression_heatmap.jpg", width=300,height=100, units = "mm", res=1000)
# # print(ht)
# # dev.off()
# 
# ctcs <- readRDS("data/cytof_objects/all_samples_ctcs_with_subtype.rds")
# assay(ctcs, "exprs") <- t(scale(t(assay(ctcs, "exprs"))))
# curr_ht <- create_expression_heatmap(ctcs, "subtype", rownames(curr_sce), "all CTCs", scale=F)
# curr_ht
