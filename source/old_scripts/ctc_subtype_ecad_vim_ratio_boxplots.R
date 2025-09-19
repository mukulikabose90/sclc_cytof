ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype_3.rds")

ctcs <- ctcs[,!is.na(ctcs$subtype)]

y <- assay(ctcs, "exprs")

df <- data.frame(t(y), colData(ctcs), check.names = FALSE)

temp <- df %>% 
  select("E-Cad","Vimentin","subtype") %>% 
  mutate("ratio" = log(`E-Cad`/Vimentin))


ggplot(temp)+
  geom_boxplot(aes(x=subtype,y=ratio,fill=subtype))



temp <- df %>% 
  select("EpCAM","Vimentin","subtype") %>% 
  mutate("ratio" = log(EpCAM/Vimentin))

ggplot(temp)+
  geom_boxplot(aes(x=subtype,y=ratio,fill=subtype))


temp <- df %>% 
  select("E-Cad","EpCAM","Vimentin","subtype") %>% 
  mutate("ratio" = log((`E-Cad`+EpCAM)/Vimentin))

ggplot(temp)+
  geom_boxplot(aes(x=subtype,y=ratio,fill=subtype))
