source("source/cytof_de_function.R")

script_seed <- 42

ctcs <- readRDS("data/cytof_objects/drug_treatment_experiment_sce_object.rds")

ctcs <- ctcs[rowData(ctcs)$marker_class == "state",]

# EMT_proteins <- c("E-Cad","SLUG","Vimentin","Twist","EpCAM","MUC-1","CD24",'CD44')
# 
# rownames(ctcs)
# 
# ctcs <- ctcs[EMT_proteins,]





y <- assay(ctcs, "exprs")


# y <- t(scale(t(y)))
dim(y)

df <- data.frame(t(y), colData(ctcs), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")
gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(ctcs)))





temp <- gg_df %>% 
  group_by(day,antigen) %>% 
  summarize(med_expr = median(expression)) %>% 
  group_by(day,antigen) %>% 
  pivot_wider(names_from = day, values_from = med_expr) %>% 
  column_to_rownames("antigen")


temp <- apply(temp, 2, FUN = function(x) x/temp[,1])


temp <- temp[,c(1,4,5,2,3)]


log_temp <- log(temp)



temp2 <- as.data.frame(log_temp) %>% 
  rownames_to_column("antigen") %>% 
  pivot_longer(!antigen, names_to = "day", values_to = "expr")


temp2$day <- factor(temp2$day, levels=c(0,3,7,10,14))

ggplot(temp2)+
  geom_line(aes(x=day,y=expr,group=antigen,color=antigen),linewidth=2)+
  facet_wrap(~antigen)+
  xlab("Day")+
  ylab("log(FC) Expression Relative to Control")




