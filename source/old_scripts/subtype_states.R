source("source/cytof_de_function.R")

script_seed <- 42

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype_3.rds")

ctcs <- ctcs[rowData(ctcs)$marker_class == "state",]

subtype_tfs <- c("POU2F3", "ASCL1","NeuroD1")

dim(ctcs)

ctcs <- ctcs[subtype_tfs,]

ctcs <- CATALYST::cluster(ctcs, features = "state",
                          xdim = 10, ydim = 10, maxK = 20, seed = script_seed)


ctcs <- runDR(ctcs, "UMAP", cells = 5e3, features = "state")


ctcs@metadata$delta_area

plotDR(ctcs, "UMAP", color_by = "meta9",scale = T)+
  geom_point(size=1)

plotDR(ctcs, "UMAP", color_by = "treatment_status",scale = T)+
  geom_point(size=1)

colData(ctcs)$new_clusters <- cluster_ids(ctcs, "meta4")
################################################################################

y <- assay(ctcs, "exprs")

df <- data.frame(t(y), colData(ctcs), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")
gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(ctcs)))




gg_df <- gg_df %>% 
  dplyr::filter(antigen %in% c("ASCL1","NeuroD1","POU2F3"))

ggplot(gg_df)+
  geom_boxplot(aes(x=new_clusters, y=expression))+
  facet_wrap(~antigen, scales="free_y")+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))




cluster_prop_df <- as.data.frame(colData(ctcs)) %>%
  group_by(new_clusters,subtype) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))%>% 
  group_by(new_clusters) %>% 
  mutate(total = sum(n))


cluster_prop_df <- cluster_prop_df %>% 
  dplyr::filter(total > 5)

ggplot(cluster_prop_df,aes(x=new_clusters, y=freq, fill=subtype))+
  geom_col()+
  geom_text(aes(y = .5,label=total))+
  xlab("Patient")+
  ylab("Percentage")+
  theme(axis.text.x = element_text(angle=70,hjust=1,vjust=1))



cluster_prop_df <- as.data.frame(colData(ctcs)) %>%
  group_by(treatment_status,new_clusters) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))%>% 
  group_by(treatment_status) %>% 
  mutate(total = sum(n))


cluster_prop_df <- cluster_prop_df %>% 
  dplyr::filter(total > 5)

ggplot(cluster_prop_df,aes(x=treatment_status, y=freq, fill=new_clusters))+
  geom_col(position = "dodge")+
  geom_text(aes(y = .5,label=total))+
  xlab("Patient")+
  ylab("Percentage")+
  theme(axis.text.x = element_text(angle=70,hjust=1,vjust=1))


#############################################

# debugonce(cytof_de)
df <- cytof_de(ctcs, method = "wilcox", metric = "mean", ident = "new_clusters")

plot_df <- df %>% 
  mutate(significance = factor(ifelse(p_adj < .05, "signif", "ns"), levels=c("signif","ns"))) %>% 
  mutate(protein = reorder_within(protein, logfc, ident_list))



# Plot DE results of significant clusters
ggplot(plot_df) +
  geom_col(aes(x=as.numeric(logfc), y=protein, fill=significance))+
  facet_wrap(~ident_list, scales = "free_y", nrow=2)+
  scale_y_reordered()+
  ylab("Protein")+
  xlab("Mean log2FC")
