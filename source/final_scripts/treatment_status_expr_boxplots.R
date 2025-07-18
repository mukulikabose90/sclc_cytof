source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam", "E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist", "SLUG", "PD-L1", "p-YAP", "CD44", "CD24")
################################################################################
# Select patients to use
################################################################################
# patients_to_use <- ctcs %>%
#   colData() %>%
#   as.data.frame() %>%
#   filter(is.na(tarla)) %>%
#   select(patient_id,treatment_status) %>%
#   distinct() %>%
#   count(patient_id) %>%
#   filter(n > 1) %>%
#   pull(patient_id) %>%
#   as.character()
# 
# # Remove samples with too few cells
# patients_to_remove <- ctcs %>%
#   colData() %>%
#   as.data.frame() %>%
#   count(patient_id,sample_num) %>%
#   filter(n < 10) %>%
#   pull(patient_id) %>%
#   unique() %>%
#   as.character()
# 
# patients_to_use <- patients_to_use[!patients_to_use %in% patients_to_remove]
# 
# curr_data <- ctcs[,ctcs$patient_id %in% patients_to_use]
################################################################################
sce <- ctcs

y <- assay(sce, "exprs")

df <- data.frame(t(y), colData(sce), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")

gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(sce)))

plot_df <- gg_df %>%
  dplyr::filter(antigen %in% markers_to_use)

plot_df$antigen <- factor(plot_df$antigen, levels = markers_to_use)

plot_df$treatment_status <- ifelse(plot_df$treatment_status == "naive", "Naive", "Treated")

plot_df$treatment_status <- factor(plot_df$treatment_status, levels = c("Naive", "Treated"))

# Remove post-tarla samples
plot_df <- plot_df %>%
  filter(tarla != "post" | is.na(tarla))

unique(plot_df$tarla)
################################################################################

stat.test <- plot_df %>%
  group_by(antigen) %>%
  wilcox_test(expression ~ treatment_status) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test


p <- ggviolin(plot_df, x="treatment_status" ,y="expression", fill="treatment_status", lwd=.3, outlier.size = .1,draw_quantiles =0.5)+
  # stat_compare_means(aes(group = treatment_status), label = "p.signif", label.x.npc = "center", label.y = 5.5,size=4.5)+
  facet_wrap(~antigen,nrow=2)+
  ylim(0,9)+
  labs(y="Expression",
       x= "")+
  scale_fill_manual(name = "Subtype",values=c("#E63946","#457B9D"))+
  theme(axis.title = element_text(size=14),
        axis.text.x = element_text(size=10,angle=45,hjust=1,vjust=1),
        strip.text = element_text(face = "bold", size=12), 
        strip.background = element_blank())+
  rremove("legend")


stat.test <- stat.test %>% add_xy_position(x = "treatment_status")
p <- p + stat_pvalue_manual(stat.test, label = "p.adj.signif")



jpeg(glue("figures/treatment_status_expression_violinplot.jpg"), width=360,height=140, units = "mm", res=1000)
print(p)
dev.off()

################################################################################
################################################################################
################################################################################
################################################################################
# 
# curr_data$treatment_status <- factor(curr_data$treatment_status, levels = c("naive", "treated"))
# 
# markers_to_use <- c("NeuroD1","ASCL1","POU2F3","p-Rb")
# 
# markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam", "E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist", "SLUG", "PD-L1", "p-YAP", "CD44", "CD24")
# p1 <- create_marker_boxplots(curr_data,markers_to_use,"patient_id",fill="treatment_status")
# 
# p1
# 
# p2 <- create_marker_boxplots(curr_data,markers_to_use,"treatment_status",fill="treatment_status")
# 
# p2
# 
# p1 <- p1+
#   rremove("legend")+
#   labs(x="Subtype",
#        y="Expression")+
#   scale_fill_manual(
#     values = c("A" = "#E57373","N" = "#FFB74D","P"= "#81C784", "I" = "#64B5F6"))
# 
# 
# 
# 
# sce <- curr_data
# 
# y <- assay(sce, "exprs")
# 
# df <- data.frame(t(y), colData(sce), check.names = FALSE)
# 
# value <- ifelse("exprs" == "exprs", "expression", "exprs")
# 
# gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
#               id.vars = names(colData(sce)))
# 
# 
# plot_df <- gg_df %>%
#   dplyr::filter(antigen %in% markers_to_use)
# 
# plot_df$antigen <- factor(plot_df$antigen, levels = markers_to_use)
# 
# p <- ggboxplot(plot_df, x="patient_id" ,y="expression", fill="treatment_status", lwd=.3, outlier.size = .1)+
#   stat_compare_means(aes(group = treatment_status), label = "p.signif",label.y = 5.5, size=3)
# 
# 
# p
# p1 <- p + 
#   facet_wrap(~antigen,nrow=3)+
#   ylim(0,6)+
#   labs(y="Expression",
#        x= "Patient ID")+
#   theme(axis.title = element_text(size=14),
#         axis.text.x = element_text(size=10,angle=45,hjust=1,vjust=1),
#         strip.text = element_text(face = "bold", size=12), 
#         strip.background = element_blank())
#   # rremove("legend")
# 
# 
# p1
# 
# jpeg(glue("figures/treatment_status_expression_patient_boxplot.jpg"), width=160,height=160, units = "mm", res=1000)
# print(p1)
# dev.off()
# 
# 
# 
# 
# p <- ggboxplot(plot_df, x="treatment_status" ,y="expression", fill="treatment_status", lwd=.3, outlier.size = .1)+
#   stat_compare_means(aes(group = treatment_status), label = "p.signif", label.x.npc = "center", label.y = 5.5,
#                      size=4.5)
# 
# 
# 
# 
# p
# p2 <- p + 
#   facet_wrap(~antigen,nrow=3)+
#   ylim(0,6)+
#   labs(y="Expression",
#        x= "Treatment Status")+
#   theme(axis.title = element_text(size=14),
#         axis.text.x = element_text(size=10,angle=45,hjust=1,vjust=1),
#         strip.text = element_text(face = "bold", size=12), 
#         strip.background = element_blank())+
#   rremove("legend")
# 
# 
# p2
# 
# jpeg(glue("figures/treatment_status_expression_boxplot.jpg"), width=160,height=160, units = "mm", res=1000)
# print(p2)
# dev.off()




