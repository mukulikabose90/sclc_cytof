source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################

markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam", "E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist", "SLUG", "PD-L1", "p-YAP", "CD44", "CD24")

col_fun = colorRamp2(c(-2, -1, 0, 1, 2), 
                     c("#313695",  # deep blue
                       "#74add1",  # light blue
                       "#f7f7f7",  # white (center, 0)
                       "#f46d43",  # light red
                       "#a50026"))

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

samples_to_use <- c("SC443-1","SC443-1P")

curr_data <- ctcs

#Scale expression
assay(curr_data, "exprs") <- t(scale(t(assay(curr_data, "exprs"))))

curr_data <- curr_data[,curr_data$collection_id %in% samples_to_use]

curr_heatmap <- matrix(NA,ncol=2,nrow=length(markers_to_use))
for(curr_sample in samples_to_use){
  
  i <- which(curr_sample == samples_to_use)
  
  y <- assay(curr_data[,curr_data$collection_id == curr_sample], "exprs")
  
  y <- y[markers_to_use,]
  
  curr_heatmap[,i] <- rowMedians(y)
  
  
}

colnames(curr_heatmap) <- samples_to_use
rownames(curr_heatmap) <- markers_to_use

curr_ht <- Heatmap(curr_heatmap, cluster_rows = F, cluster_columns = F,column_names_rot = 0,
                   name="Median Scaled\n   Expression", column_title = "", col = col_fun)


p1 <- draw(curr_ht)
p1

jpeg(glue("figures/pleural_effusion_expression_heatmap.jpg"), width=260,height=150, units = "mm", res=1000)
print(p1)
dev.off()
