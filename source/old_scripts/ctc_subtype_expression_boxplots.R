source("source/cytof_de_function.R")

script_seed <- 42
set.seed(script_seed)
################################################################################


ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype_3.rds")

y <- assay(ctcs, "exprs")

dim(y)

# y <- t(scale(t(y)))

# Scale data. Scale protein expression between 0-1 across all cells
# scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
# y <- t(apply(y, MARGIN = 1, FUN = function(X) scale_values(X)))

df <- data.frame(t(y), colData(ctcs), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")
gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(ctcs)))


gg_df <- gg_df %>% 
  dplyr::filter(!is.na(subtype))

# ggplot(gg_df)+
#   geom_boxplot(aes(x=subtype, y=expression))+
#   facet_wrap(~antigen, scales="free_y")+
#   theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))


p <- ggboxplot(gg_df, x = "subtype", y = "expression",
               color = "subtype", palette = "jco", facet.by = "antigen",
               scales="free_y")


p + stat_compare_means(method = "anova")+
  ylab("Scaled Expression")




