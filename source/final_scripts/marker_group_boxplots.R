
sce <- readRDS("data/cytof_objects/sclc_all_samples_with_clusters.rds")

markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam", "E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist", "SLUG", "PD-L1", "p-YAP", "CD44", "CD24")

marker_group <- list()
marker_group[["SCLC TFs"]] <- c("ASCL1","POU2F3","NeuroD1")
marker_group[["Therapeutic targets"]] <- c("DLL3","CD24","PD-L1")
marker_group[["EMT Markers"]] <- c("E-Cad","Vimentin","EpCAM","MUC-1","CD44","SLUG","Twist","Alcam")
marker_group[["Resistance"]] <- c("p-YAP")



marker_groups <- c()
for(i in markers_to_use){
  marker_groups <- append(marker_groups, names(marker_group)[sapply(marker_group, FUN = function(x) i %in% x)])
}

names(markers_to_use) <- marker_groups

marker_df <- markers_to_use %>% 
  enframe()

colnames(marker_df) <- c("marker_group","antigen")



y <- assay(sce, "exprs")

df <- data.frame(t(y), colData(sce), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")

gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(sce)))




################################################################################

plot_df <- gg_df %>%
  dplyr::filter(antigen %in% markers_to_use)

dim(plot_df)


plot_df <- merge(plot_df,marker_df,by="antigen",all=T)

plot_df$antigen <- factor(plot_df$antigen, levels= markers_to_use)

p <- ggboxplot(plot_df, x="new_clusters" ,y="expression")

p + 
  facet_wrap(~antigen, scales = "free")+
  stat_compare_means()


colnames(plot_df)

plots <- list()
for(curr_group in unique(marker_groups)){
  curr_df <- plot_df %>% 
    filter(marker_group == curr_group)
  
  p <- ggboxplot(curr_df, x="new_clusters" ,y="expression")
  
  p1 <- p + 
    facet_grid(marker_group~antigen, scales = "free")+
    stat_compare_means()
  
  plots <- append(plots,list(p1))
  
}


ggarrange(plotlist = plots,ncol=1)


