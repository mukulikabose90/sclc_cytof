source("source/cytof_de_function.R")

sce <- readRDS("data/cytof_objects/sclc_cytof_sce_object_with_clusters.rds")

sce <- sce[,colData(sce)$condition == "cancer"]

ctcs <- sce[,colData(sce)$new_clusters %in% c(1,2)]
# ctcs <- sce

ctcs <- ctcs[rowData(ctcs)$marker_class == "state",]


y <- assay(ctcs, "exprs")


dim(y)

# y <- y[rownames(y) %in% c("POU2F3","NeuroD1","ASCL1"),]

cor_mat <- cor(t(y), method = "spearman")



col_fun2 = colorRamp2(c(-.5, 0, .5), c("blue", "white", "red"))
Heatmap(cor_mat, col = col_fun2, name='Spearman\nCorrelation')



plot(t(y)[,1],t(y)[,1])

plot(t(y)[,1],t(y)[,2])
