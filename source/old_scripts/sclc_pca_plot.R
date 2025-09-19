

sce <- readRDS("data/cytof_objects/sclc_cytof_sce_object.rds")

sce <- sce[,colData(sce)$condition != "cell_line"]

pbMDS(sce, color_by = "condition", label_by = "patient_id")


y <- assay(sce, "exprs")




pca_mat <- prcomp(t(y), scale = FALSE)



temp <- data.frame(pca_mat$x) %>% 
  select(PC1, PC2) %>% 
  cbind(.,colData(sce)) %>% 
  group_by(patient_id,condition) %>% 
  summarize(mean(PC1),mean(PC2))

colnames(temp)[3:4] <- c("PC1","PC2")

ggplot(temp)+
  geom_point(aes(x=PC1,y=PC2,color=condition),size=3)
