source("source/cytof_de_function.R")

script_seed <- 42

ctcs <- readRDS("data/cytof_objects/drug_treatment_experiment_sce_object.rds")

ctcs <- ctcs[rowData(ctcs)$marker_class == "state",]

EMT_proteins <- c("E-Cad","SLUG","Vimentin","Twist","EpCAM","MUC-1","CD24",'CD44')

rownames(ctcs)

ctcs <- ctcs[EMT_proteins,]

ctcs <- CATALYST::cluster(ctcs, features = "state",
                          xdim = 10, ydim = 10, maxK = 20, seed = script_seed)


ctcs <- runDR(ctcs, "UMAP", cells = 5e3, features = "state")


ctcs@metadata$delta_area

plotDR(ctcs, "UMAP", color_by = "meta10",scale = T)+
  geom_point(size=1)


colData(ctcs)$new_clusters <- cluster_ids(ctcs, "meta10")

################################################################################

################################################################################

cluster_prop_df <- as.data.frame(colData(ctcs)) %>%
  group_by(new_clusters,day) %>%
  summarise(n = n()) %>%
  group_by(day) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n / total*100))

cluster_prop_df <- cluster_prop_df %>% 
  dplyr::filter(total > 5)


temp <- cluster_prop_df %>% 
  select(day,new_clusters,freq) %>% 
  group_by(day,new_clusters) %>% 
  pivot_wider(names_from = day, values_from = freq) %>% 
  column_to_rownames("new_clusters")


temp <- apply(temp, 2, FUN = function(x) x/temp[,1])

temp <- temp[,c(1,4,5,2,3)]


# log_temp <- log(temp)


temp2 <- as.data.frame(temp) %>% 
  rownames_to_column("new_clusters") %>% 
  pivot_longer(!new_clusters, names_to = "day", values_to = "freq")


temp2$day <- factor(temp2$day, levels=c(0,3,7,10,14))

ggplot(temp2)+
  geom_line(aes(x=day,y=freq,group=new_clusters,color=new_clusters),linewidth=2)+
  facet_wrap(~new_clusters)+
  xlab("Day")+
  ylab("log(FC) Expression Relative to Control")













