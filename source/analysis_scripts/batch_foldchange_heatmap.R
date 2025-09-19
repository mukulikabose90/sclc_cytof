
source("source/sclc_cytof_functions.R")

sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  filter(collection_id == "SC414-1") %>% 
  pull(experiment_id) %>% 
  unique()


sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  filter(collection_id == "NORMAL1-1") %>% 
  pull(experiment_id) %>% 
  unique()

temp <- sce[,sce$sample_id == "SC414-1_515600"|sce$sample_id == "NORMAL1-1_515600"|sce$sample_id == "SC414-1_513549"|sce$sample_id == "NORMAL1-1_513549"]

# sce[,sce$sample_id == "NORMAL1-1_515600"]


temp@assays@data$exprs

y <- assay(temp, "exprs")

df <- data.frame(t(y), colData(temp), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")

gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(temp)))


markers_to_use <- readRDS("data/state_markers.rds")
plot_df <- gg_df %>%
  dplyr::filter(antigen %in% markers_to_use)


ggplot(plot_df)+
  geom_density(aes(x=expression,color=experiment_id))+
  facet_wrap(~antigen, scales="free")+
  theme_classic()


sce <- readRDS("data/cytof_objects/sclc_all_samples_object_no_qc.rds")

temp1 <- sce[,sce$sample_id == "SC414-2_515600"]
temp2 <- sce[,sce$sample_id == "NORMAL1-1_515600"]
temp3 <- sce[,sce$sample_id == "SC414-2_513549"]
temp4 <- sce[,sce$sample_id == "NORMAL1-1_513549"]

col1 <- rowMeans(temp2@assays@data$exprs[markers_to_use,])/rowMeans(temp2@assays@data$exprs[markers_to_use,])
col2 <- rowMeans(temp1@assays@data$exprs[markers_to_use,])/rowMeans(temp2@assays@data$exprs[markers_to_use,])
col3 <- rowMeans(temp4@assays@data$exprs[markers_to_use,])/rowMeans(temp4@assays@data$exprs[markers_to_use,])
col4 <- rowMeans(temp3@assays@data$exprs[markers_to_use,])/rowMeans(temp4@assays@data$exprs[markers_to_use,])



df1 <- data.frame(names(col1),col1,"10/10 ctrl")
colnames(df1) <- c("marker","fc","date")
rownames(df1) <- NULL

df2 <- data.frame(names(col2),col2,"10/10")
colnames(df2) <- c("marker","fc","date")
rownames(df2) <- NULL

df3 <- data.frame(names(col3),col3,"9/12 ctrl")
colnames(df3) <- c("marker","fc","date")
rownames(df3) <- NULL

df4 <- data.frame(names(col4),col4,"9/12")
colnames(df4) <- c("marker","fc","date")
rownames(df4) <- NULL


df <- rbind(df1,df2,df3,df4)

df$marker

ggplot(df,aes(x=marker,y=fc,color=date))+
  geom_point()+
  geom_line(aes(group=date))


heatmap <- cbind(col1,col2,col3,col4)
colnames(heatmap) <- c("10/10 ctrl", "10/10","9/12 ctrl", "9/12")

col_fun = colorRamp2(c(-2, -1, 0, 1, 2), c("white", "blue", "black", "yellow", "white"))
Heatmap(log2(heatmap),cluster_columns = F,cluster_rows=F, name="log2(FC)",col=col_fun)

dev.off()
col1 <- rowMedians(temp2@assays@data$exprs[markers_to_use,])/rowMedians(temp2@assays@data$exprs[markers_to_use,])
col2 <- rowMedians(temp1@assays@data$exprs[markers_to_use,])/rowMedians(temp2@assays@data$exprs[markers_to_use,])
col3 <- rowMedians(temp4@assays@data$exprs[markers_to_use,])/rowMedians(temp4@assays@data$exprs[markers_to_use,])
col4 <- rowMedians(temp3@assays@data$exprs[markers_to_use,])/rowMedians(temp4@assays@data$exprs[markers_to_use,])
heatmap <- cbind(col1,col2,col3,col4)
colnames(heatmap) <- c("10/10 ctrl", "10/10","9/12 ctrl", "9/12")

rownames(heatmap) <- markers_to_use

col_fun = colorRamp2(c(-2, -1, 0, 1, 2), c("white", "blue", "black", "yellow", "white"))
Heatmap(log2(heatmap),cluster_columns = F,cluster_rows=F, name="log2(FC)",col=col_fun)

