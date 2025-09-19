library(reshape2)
library(ggpubr)
# This script compares the protein expression of the same sample across different 
# CyTOF experiments

sce <- readRDS("data/cytof_objects/sclc_cytof_sce_object.rds")

sce <- sce[,colData(sce)$sample_type != "cell_line"]

# Add cell IDs
colData(sce)$cell_id <- 1:ncol(sce)
y <- assay(sce, "exprs")
df <- data.frame(t(y), colData(sce), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")
gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(sce))) 

my_comparisons <- list( c("513549", "514521"), c("513549", "515600"), c("515600", "514521") )

curr_patient <- 414

ggboxplot(gg_df, x = "experiment_id", y = "expression",
          color = "condition", palette = "jco",
          add = "jitter", facet.by = c("patient_id","antigen")) + 
  
  stat_compare_means(method = "wilcox.test", comparisons = my_comparisons)



compare_means(data = gg_df, x = "experiment_id", y = "expression",
                   facet.by = c("patient_id","antigen"),
method = "wilcox.test", comparisons = my_comparisons)
