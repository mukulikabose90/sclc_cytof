source("source/cytof_de_function.R")

script_seed <- 42
set.seed(script_seed)

metric_to_use <- "median"
################################################################################
# Differential expression

# Set significantly cancer abundant clusters
signif_clusters <- c(3,4)

sce <- readRDS("data/cytof_objects/sclc_first_draw_with_clusters.rds")

df <- cytof_de(sce, method = "wilcox", metric = metric_to_use, ident = "new_clusters")


plot_df <- df %>% 
  mutate(significance = factor(ifelse(p_adj < .05, "*", ""), levels=c("*",""))) %>% 
  mutate(protein = reorder_within(protein, logfc, ident_list)) %>% 
  mutate(star_x = ifelse(logfc > 0, logfc+.25,logfc-.5)) %>% 
  mutate(up_down= ifelse(logfc>0,"up","down"))


plot_df$ident_list <- paste0("Cluster ",plot_df$ident_list)
plot_df$ident_list <- factor(plot_df$ident_list, levels=paste0("Cluster ", c(1:10)))


signif_df <- plot_df %>% 
  dplyr::filter(ident_list %in% paste0("Cluster ", signif_clusters))


x_axis_label <- gsub("m","M",metric_to_use)

p1 <- ggplot(signif_df,aes(x=as.numeric(logfc), y=protein, fill=up_down))+
  geom_col(color="darkgray",size=.001)+
  geom_text(aes(x=star_x, label=significance), size = 5)+
  facet_wrap(~ident_list, scales = "free_y",nrow=1)+
  scale_y_reordered()+
  labs(fill = "")+
  guides(fill="none")+
  scale_fill_manual(values=c("royalblue", "firebrick1"))+
  ylab("Protein")+
  xlab(glue("{x_axis_label} log(FC)"))+
  theme_classic()+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=20), 
        axis.text = element_text(color = "black", size=15),
        axis.title = element_text(size=20),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18))


p2 <- ggplot(plot_df,aes(x=as.numeric(logfc), y=protein, fill=up_down))+
  geom_col(color="darkgray",size=.001)+
  geom_text(aes(x=star_x, label=significance), size = 5)+
  facet_wrap(~ident_list, scales = "free_y", nrow=2)+
  scale_y_reordered()+
  labs(fill = "")+
  guides(fill="none")+
  scale_fill_manual(values=c("royalblue", "firebrick1"))+
  ylab("Protein")+
  xlab(glue("{x_axis_label} log(FC)"))+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=20), 
        axis.text = element_text(color = "black", size=15),
        axis.title = element_text(size=20),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18))

p1


jpeg(glue("figures/first_draw_cluster_{metric_to_use}_de_signif.jpg"), width=120,height=100, units = "mm", res=1000)
print(p1)
dev.off()

jpeg(glue("figures/first_draw_cluster_{metric_to_use}_de_all.jpg"), width=120,height=100, units = "mm", res=1000)
print(p2)
dev.off()









