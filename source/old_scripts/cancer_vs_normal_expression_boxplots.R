source("source/cytof_de_function.R")


sce <- readRDS("data/cytof_objects/sclc_cytof_sce_object.rds")

sce <- sce[,colData(sce)$condition != "cell_line"]


y <- assay(sce, "exprs")

# Scale data. Scale protein expression between 0-1 across all cells
# scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
# y <- t(apply(y, MARGIN = 1, FUN = function(X) scale_values(X)))

df <- data.frame(t(y), colData(sce), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")
gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(sce)))

ggplot(gg_df)+
  geom_boxplot(aes(x=condition, y=expression))+
  facet_wrap(~antigen, scales="free_y")+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))


p <- ggboxplot(gg_df, x = "condition", y = "expression",
               color = "condition", palette = "jco", facet.by = "antigen")


p + stat_compare_means(method = "wilcox.test",label.y = 7.5)
