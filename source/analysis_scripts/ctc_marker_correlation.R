

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

y <- assay(ctcs, "exprs")

df <- data.frame(t(y), colData(ctcs), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")

gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(ctcs)))

########################################################################
markers_to_use <- c("DLL3","ASCL1")

plot_df <- gg_df %>% 
  filter(antigen %in% markers_to_use) %>% 
  select(collection_id,antigen,expression) %>% 
  group_by(collection_id, antigen) %>% 
  mutate(med_expr = median(expression)) %>%
  select(collection_id,antigen,med_expr) %>% 
  distinct() %>% 
  pivot_wider(names_from = antigen, values_from = med_expr)
  

ggscatter(plot_df,x ="ASCL1",y="DLL3",add = "reg.line")+
  stat_cor(label.y = 3) 


plot_df <- gg_df %>% 
  filter(antigen %in% markers_to_use) %>% 
  select(cell_id,antigen,expression) %>%
  distinct() %>% 
  pivot_wider(names_from = antigen, values_from = expression)

ggscatter(plot_df,x ="ASCL1",y="DLL3",add = "reg.line")+
  stat_cor(label.y = 7) 

########################################################################
markers_to_use <- c("POU2F3","ASCL1")

plot_df <- gg_df %>% 
  filter(antigen %in% markers_to_use) %>% 
  select(collection_id,antigen,expression) %>% 
  group_by(collection_id, antigen) %>% 
  mutate(med_expr = median(expression)) %>%
  select(collection_id,antigen,med_expr) %>% 
  distinct() %>% 
  pivot_wider(names_from = antigen, values_from = med_expr)


ggscatter(plot_df,x ="ASCL1",y="POU2F3",add = "reg.line")+
  stat_cor(label.y = 3) 

plot_df <- gg_df %>% 
  filter(antigen %in% markers_to_use) %>% 
  select(cell_id,antigen,expression,treatment_status) %>%
  distinct() %>% 
  pivot_wider(names_from = antigen, values_from = expression)

ggscatter(plot_df,x ="ASCL1",y="POU2F3",add = "reg.line")+
  facet_wrap(~treatment_status)+
  stat_cor(label.y = 7) 

########################################################################
markers_to_use <- c("NeuroD1","ASCL1")

plot_df <- gg_df %>% 
  filter(antigen %in% markers_to_use & expression > 1) %>% 
  select(collection_id,antigen,expression) %>% 
  group_by(collection_id, antigen) %>% 
  mutate(med_expr = median(expression)) %>%
  select(collection_id,antigen,med_expr) %>% 
  distinct() %>% 
  pivot_wider(names_from = antigen, values_from = med_expr)


ggscatter(plot_df,x ="ASCL1",y="NeuroD1",add = "reg.line")+
  stat_cor(label.y = 3) 

plot_df <- gg_df %>% 
  filter(antigen %in% markers_to_use & expression > 1) %>% 
  select(cell_id,antigen,expression,treatment_status) %>%
  pivot_wider(names_from = antigen, values_from = expression)

ggscatter(plot_df,x ="ASCL1",y="NeuroD1",add = "reg.line")+
  facet_wrap(~treatment_status)+
  stat_cor(label.y = 7) 

