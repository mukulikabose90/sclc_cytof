source("source/cytof_de_function.R")

script_seed <- 42
set.seed(script_seed)
################################################################################
# Differential expression

metric_to_use <- "median"
ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

length(as.character(unique(ctcs$patient_id)))
dim(ctcs)

ctcs <- ctcs[,!is.na(ctcs$subtype)]

df <- cytof_de(ctcs, method = "wilcox", metric = metric_to_use, ident = "subtype")

plot_df <- df %>% 
  mutate(significance = factor(ifelse(p_adj < .05, "*", ""), levels=c("*",""))) %>% 
  mutate(protein = reorder_within(protein, logfc, ident_list)) %>% 
  mutate(star_x = ifelse(logfc > 0, logfc+.25,logfc-.5)) %>% 
  mutate(up_down= ifelse(logfc>0,"up","down"))

plot_df$ident_list <- ifelse(plot_df$ident_list == "A", "ASCL1",ifelse(plot_df$ident_list == "P", "POU2F3",ifelse(plot_df$ident_list == "N", "NeuroD1","Inflamed")))

x_axis_label <- gsub("m","M",metric_to_use)

p1 <- ggplot(plot_df,aes(x=as.numeric(logfc), y=protein, fill=up_down))+
  geom_col(color="darkgray",size=.001)+
  geom_text(aes(x=star_x, label=significance), size = 5)+
  facet_wrap(~ident_list, scales = "free_y")+
  scale_y_reordered()+
  labs(fill = "")+
  guides(fill="none")+
  scale_fill_manual(values=c("royalblue", "firebrick1"))+
  ylab("Protein")+
  xlab(glue("{x_axis_label} log(FC)"))+
  theme_classic()+
  theme(panel.grid.minor = element_blank(), 
        plot.title = element_text(size=30),
        strip.text = element_text(face = "bold", size=20), 
        axis.text = element_text(color = "black", size=15),
        axis.title = element_text(size=20),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18))


png("figures/ctcs_subtype_de.png", width = 12, height = 14, units = "in", res = 1200)
# tiff("figures/ctcs_subtype_heatmap.tiff", width = 12, height = 5, units = "in", res = 1200)
print(p1)
# dev.off()
dev.off()

