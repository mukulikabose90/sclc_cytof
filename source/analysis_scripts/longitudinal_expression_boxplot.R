source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")


patients_to_use <- ctcs %>% 
  colData() %>% 
  as.data.frame() %>% 
  select(patient_id,collection_id) %>% 
  distinct() %>% 
  count(patient_id) %>% 
  filter(n > 1) %>% 
  pull(patient_id)

curr_ctcs <- ctcs[,ctcs$patient_id %in% patients_to_use]

markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam", "E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist", "SLUG", "PD-L1", "p-YAP", "CD44", "CD24")

sce <- curr_ctcs

y <- assay(sce, "exprs")

df <- data.frame(t(y), colData(sce), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")

gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(sce)))


plot_df <- gg_df %>%
  dplyr::filter(antigen %in% markers_to_use)

plot_df$antigen <- factor(plot_df$antigen, levels = markers_to_use)

plot_df <- plot_df %>% 
  group_by(patient_id) %>% 
  mutate(sample_new = as.integer(factor(sample_num))) 

class(plot_df$sample_new)

plot_df$sample_new <- as.character(plot_df$sample_new)

p <- ggboxplot(plot_df, x="patient_id" ,y="expression", fill="sample_new", lwd=.3, outlier.size = .1)+
  stat_compare_means(aes(group = sample_new), label = "p.signif",label.y = 5.5, size=3)

p
p1 <- p + 
  facet_wrap(~antigen,nrow=3)+
  ylim(0,6)+
  labs(y="Expression",
       x= "Patient ID")+
  theme(axis.title = element_text(size=14),
        axis.text.x = element_text(size=10,angle=45,hjust=1,vjust=1),
        strip.text = element_text(face = "bold", size=12), 
        strip.background = element_blank())+
  rremove("legend")


p1


jpeg(glue("figures/longitudinal_expression_boxplot.jpg"), width=350,height=150, units = "mm", res=1000)
print(p1)
dev.off()

################################################################################

patients_to_use <- ctcs %>% 
  colData() %>% 
  as.data.frame() %>% 
  filter(is.na(tarla)) %>% 
  select(patient_id,collection_id) %>% 
  distinct() %>% 
  count(patient_id) %>% 
  filter(n > 1) %>% 
  pull(patient_id)

curr_ctcs <- ctcs[,ctcs$patient_id %in% patients_to_use]

markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam", "E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist", "SLUG", "PD-L1", "p-YAP", "CD44", "CD24")

sce <- curr_ctcs

y <- assay(sce, "exprs")

df <- data.frame(t(y), colData(sce), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")

gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(sce)))


plot_df <- gg_df %>%
  dplyr::filter(antigen %in% markers_to_use)

plot_df$antigen <- factor(plot_df$antigen, levels = markers_to_use)

plot_df <- plot_df %>% 
  group_by(patient_id) %>% 
  mutate(sample_new = as.integer(factor(sample_num))) 

class(plot_df$sample_new)

plot_df$sample_new <- as.character(plot_df$sample_new)

p <- ggboxplot(plot_df, x="patient_id" ,y="expression", fill="sample_new", lwd=.3, outlier.size = .1)+
  stat_compare_means(aes(group = sample_new), label = "p.signif",label.y = 5.5, size=3)

p
p1 <- p + 
  facet_wrap(~antigen,nrow=3)+
  ylim(0,6)+
  labs(y="Expression",
       x= "Patient ID")+
  theme(axis.title = element_text(size=14),
        axis.text.x = element_text(size=10,angle=45,hjust=1,vjust=1),
        strip.text = element_text(face = "bold", size=12), 
        strip.background = element_blank())+
  rremove("legend")


p1



jpeg(glue("figures/longitudinal_notarla_expression_boxplot.jpg"), width=300,height=150, units = "mm", res=1000)
print(p1)
dev.off()

