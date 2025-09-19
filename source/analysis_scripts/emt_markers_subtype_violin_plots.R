source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

y <- assay(ctcs, "exprs")

df <- data.frame(t(y), colData(ctcs), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")

gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(ctcs)))

markers_to_use <- c("E-Cad","Vimentin","EpCAM","Twist","SLUG")

plot_df <- gg_df %>%
  dplyr::filter(antigen %in% markers_to_use)


plot_df$subtype <- factor(plot_df$subtype, levels = c("A","N","P","I"))



curr_plot_df <- plot_df %>% 
  filter(antigen %in% c("SLUG","Twist","Vimentin")) 

stat.test <- compare_means(group.by = "antigen",
  expression ~ subtype, data = curr_plot_df,
  method = "kruskal.test")

curr_plot_df <- stat.test %>% 
  select(antigen,p.format) %>% 
  merge(curr_plot_df, by="antigen") %>% 
  mutate(pheight = ifelse(antigen == "SLUG",2,ifelse(antigen=="Twist",4,7.5))) %>% 
  mutate(p.format = paste0("p: ",p.format))

p1 <- ggviolin(curr_plot_df, x="subtype" ,y="expression", fill="subtype", draw_quantiles = 0.5)+
  facet_wrap(~antigen,scales="free_y",ncol=3)+
  geom_text(aes(label=p.format,y=pheight),x=1.5,size=2)+
  ylim(0, NA)+
  theme_classic()+
  theme(strip.text = element_text(face = "bold", size=12), 
        strip.background = element_blank())+
  scale_fill_manual(
    values = c("A" = "#dd4b33","N" = "#F1FAEE","P"= "#A8DADC", "I" = "#457B9D"))+
  labs(y="Expression",
       x=NULL,
       fill="Subtype")+
  rremove("legend")

p1


curr_plot_df <- plot_df %>% 
  filter(!antigen %in% c("SLUG","Twist","Vimentin")) 

stat.test <- compare_means(group.by = "antigen",
                           expression ~ subtype, data = curr_plot_df,
                           method = "kruskal.test")

curr_plot_df <- stat.test %>% 
  select(antigen,p.format) %>% 
  merge(curr_plot_df, by="antigen") %>% 
  mutate(pheight = ifelse(antigen == "E-Cad",5,7)) %>% 
  mutate(p.format = paste0("p: ",p.format))

p2 <- ggviolin(curr_plot_df, x="subtype" ,y="expression", fill="subtype", draw_quantiles = 0.5)+
  facet_wrap(~antigen,scales="free_y",ncol=3)+
  geom_text(aes(label=p.format,y=pheight),x=1.5,size=2)+
  ylim(0, NA)+
  theme_classic()+
  theme(strip.text = element_text(face = "bold", size=12), 
        strip.background = element_blank())+
  scale_fill_manual(
    values = c("A" = "#dd4b33","N" = "#F1FAEE","P"= "#A8DADC", "I" = "#457B9D"))+
  labs(y="Expression",
       x="",
       fill="Subtype")

p2


jpeg("figures/subtype_emt_violin_plots.jpg", width=170,height=120, units = "mm", res=1000)
grid.arrange(grobs = lapply(
  list(p1, p2),
  set_panel_size,
  width = unit(4.5, "cm"),
  height = unit(4.5, "cm")
))
dev.off()


