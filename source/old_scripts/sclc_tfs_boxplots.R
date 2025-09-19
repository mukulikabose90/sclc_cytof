
# This script plots the protein expression distribution for the 4 SCLC subtype
# transcription factors.

sce <- readRDS("data/cytof_objects/sclc_cytof_sce_object.rds")

sce <- sce[,colData(sce)$condition != "cell_line"]

y <- assay(sce, "exprs")
df <- data.frame(t(y), colData(sce), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")
gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(sce)))


sclc_tfs <- c("NeuroD1","p-YAP","ASCL1","POU2F3")

temp <- gg_df %>% 
  dplyr::filter(antigen %in% sclc_tfs)


ggplot(temp)+
  geom_boxplot(aes(x=condition,y=expression,fill=condition))+
  facet_grid(~antigen,scales="free_y")


p <- ggboxplot(temp, x = "condition", y = "expression",
               color = "condition", palette = "jco",
               facet.by = "antigen", short.panel.labs = T,nrow = 1,
               scales = "free_y", legend.side="right")

# Use only p.format as label. Remove method name.
p + stat_compare_means(label = "p.format")+
  theme(legend.position = "right")


sum(y==0)
dim(y)
summary(t(y))
