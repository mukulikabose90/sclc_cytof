source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################
# Differential expression

metric_to_use <- "mean"
ctcs <- readRDS("data/cytof_objects/all_samples_ctcs_with_subtype.rds")

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

p1 <- ggplot(plot_df,aes(x=as.numeric(logfc), y=protein, fill=logfc))+
  geom_col(color="darkgray",size=.001)+
  geom_text(aes(x=star_x, label=significance), size = 3)+
  facet_wrap(~ident_list, scales = "free", nrow=2)+
  scale_y_reordered()+
  labs(fill = "")+
  guides(fill="none")+
  # scale_fill_manual(values=c("royalblue", "firebrick1"))+
  scale_fill_gradient2(low = "royalblue4", mid = "white", high = "firebrick", midpoint = 0)+
  ylab("Protein")+
  xlab(glue("{x_axis_label} log(FC)"))+
  theme_classic()+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=8), 
        axis.text = element_text(color = "black", size=6),
        axis.title = element_text(size=8),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8))

p1

jpeg(glue("figures/subtype_{metric_to_use}_de.jpg"), width=140,height=100, units = "mm", res=1000)
print(p1)
dev.off()

