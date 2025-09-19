source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################
################################################################################
# Get state markers
################################################################################
state_markers <- readRDS("data/state_markers.rds")
markers_to_use <- state_markers

sce <- readRDS("data/cytof_objects/sclc_all_samples_object.rds")

y <- assay(sce, "exprs")

df <- data.frame(t(y), colData(sce), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")

gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(sce)))


plot_df <- gg_df %>%
  dplyr::filter(antigen %in% markers_to_use) %>% 
  group_by(collection_id,antigen) %>% 
  mutate(median_expr = median(expression)) %>%
  mutate(mean_expr = mean(expression)) %>% 
  select(mean_expr,median_expr,collection_id,antigen,condition) %>% 
  distinct()



# p <- 
  ggviolin(plot_df, x="condition" ,y="median_expr", fill="condition", lwd=.3, outlier.size = .1)+
  stat_compare_means()+ 
  facet_wrap(~antigen,scales="free_y")+
  ylim(0,6)+
  labs(y="Median Expression",
       x= "Condition")+
  theme(axis.title = element_text(size=14),
        axis.text.x = element_text(size=10,angle=45,hjust=1,vjust=1),
        strip.text = element_text(face = "bold", size=12), 
        strip.background = element_blank())+
  rremove("legend")

  
  ggviolin(plot_df, x="condition" ,y="mean_expr", fill="condition", lwd=.3, outlier.size = .1)+
    stat_compare_means()+ 
    facet_wrap(~antigen,scales="free_y")+
    ylim(0,6)+
    labs(y="Mean Expression",
         x= "Condition")+
    theme(axis.title = element_text(size=14),
          axis.text.x = element_text(size=10,angle=45,hjust=1,vjust=1),
          strip.text = element_text(face = "bold", size=12), 
          strip.background = element_blank())+
    rremove("legend")






