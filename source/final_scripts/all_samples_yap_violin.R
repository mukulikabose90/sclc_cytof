source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################
cluster_colors <- c("#dd4b33", "#F1FAEE", "#A8DADC", "#457B9D")

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

markers_to_use <- c("p-YAP")

samples_to_use <- ctcs@colData %>% 
  as.data.frame()  %>% 
  count(collection_id,subtype,treatment_status) %>% 
  group_by(collection_id) %>% 
  mutate(total = sum(n)) %>% 
  filter(total > 10) %>% 
  pull(collection_id) %>% 
  unique() %>% 
  as.character() %>% 
  sort()

########################

curr_data <- ctcs[,ctcs$collection_id %in% samples_to_use]

y <- assay(curr_data, "exprs")

df <- data.frame(t(y), colData(curr_data), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")

gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(curr_data)))

plot_df <- gg_df %>%
  dplyr::filter(antigen %in% markers_to_use)

plot_df$antigen <- factor(plot_df$antigen, levels = markers_to_use)


# plot_df$sample_num <- ifelse(plot_df$sample_num == "4" & plot_df$patient_id == "SC338", "1", plot_df$sample_num)
# plot_df$sample_num <- ifelse(plot_df$sample_num == "1", "1", "2")

################################################################################

# stat.test <- plot_df %>%
#   group_by(patient_id) %>%
#   wilcox_test(expression ~ sample_num) %>%
#   adjust_pvalue(method = "bonferroni") %>%
#   add_significance()

length(samples_to_use)

p <- ggviolin(plot_df, x = "collection_id",y="expression", lwd=.3,
              outlier.size = .1, draw_quantiles = .5,fill="#457B9D")+
  facet_grid(. ~ treatment_status, scales = "free", space = "free") +
  ylim(0,NA)+
  labs(y="Expression",
       x= "")+
  # scale_fill_manual(name = "Sample Number",values=c("#E63946","#457B9D"))+
  theme(axis.title = element_text(size=20),
        axis.text.x = element_text(size=18,angle=45,hjust=1,vjust=1),
        strip.text = element_text(face = "bold", size=16), 
        strip.background = element_blank(),
        legend.position = "right")
rremove("legend")

p

stat.test <- stat.test %>% add_xy_position(x="patient_id")

p <- p + stat_pvalue_manual(stat.test, label = "p.adj.signif", size=5, tip.length = 0,
                            y.position = c(4,2,2.2))

tiff(glue("figures/pyap_all_samples_expression_violin.tiff"), width=260,height=150, units = "mm", res=600)
print(p)
dev.off()
  

plot_df$treatment_status <- ifelse(plot_df$treatment_status == "naive", "Naive", "Treated")

p <- ggboxplot(plot_df, x = "collection_id",y="expression", lwd=.3,
              outlier.size = 1, fill="#457B9D")+
  facet_grid(. ~ treatment_status, scales = "free", space = "free") +
  ylim(0,NA)+
  labs(y="Expression",
       x= "")+
  theme_minimal()+
  theme(axis.title = element_text(size=20),
        axis.text.x = element_text(size=18,angle=45,hjust=1,vjust=1),
        axis.text.y = element_text(size=18),
        strip.text = element_text(face = "bold", size=16), 
        strip.background = element_blank(),
        legend.position = "right",panel.border = element_rect(color = "black", fill = NA, linewidth = 1))


p
