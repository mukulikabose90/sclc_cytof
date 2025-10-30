source("source/sclc_cytof_functions.R")

set.seed(42)
################################################################################
# Read in data
################################################################################
ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

ctcs <- ctcs[,ctcs$collection_id != "MDA-SC399-2"]

cluster_colors <- c("#dd4b33", "#D1DACF", "#A8DADC", "#457B9D")

################################################################################
# Sample n cells from each patient
################################################################################
n_cells <- 30

sampled_data <- ctcs@colData %>% 
  as.data.frame() %>%
  group_by(collection_id) %>%
  filter(n() >= n_cells) %>%        # keep only patients with ≥ n_cells cells
  slice_sample(n = n_cells) %>%
  ungroup()

sampled_data$treatment_status <- ifelse(sampled_data$treatment_status == "naive","Naive","CTX ± ICI")

sampled_data$treatment_status <- ifelse(is.na(sampled_data$tarla), sampled_data$treatment_status,
                                   ifelse(sampled_data$tarla == "pre", sampled_data$treatment_status, "Tarla"))

################################################################################
# Create plot dataframe
################################################################################
plot_df <- sampled_data %>% 
  select(treatment_status,subtype,tarla) %>% 
  cbind(1) %>% 
  rename("n" = "1")

plot_df$treatment_status <- factor(plot_df$treatment_status, levels=c("Naive","CTX ± ICI","Tarla"))

plot_df$subtype <- factor(plot_df$subtype, levels=c("A","N","P","Mes"))

plot_df %>% 
  count(subtype, treatment_status) %>% 
  group_by(treatment_status) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n/total)*100) %>% 
  arrange(treatment_status) %>% 
  select(treatment_status,freq,subtype)

plot_df_long <- to_lodes_form(data.frame(plot_df),
                              key = "category", value = "group", id = "cohort",
                              axes = 1:2) %>% 
  add_count(group) %>% 
  group_by(category) %>% 
  mutate(total=sum(n)) %>% 
  mutate(freq=sprintf("%.1f",(nn/total)*100)) 


plot_df_long_left <- subset(plot_df_long, group %in% c("Naive","CTX ± ICI","Tarla"))
plot_df_long_right <- subset(plot_df_long, !group %in% c("Naive","CTX ± ICI","Tarla"))

################################################################################
# Plot alluvial 
################################################################################
p1 <- ggplot(data = plot_df_long,
             aes(x = category, stratum = group, alluvium = cohort, y = total)) +
  coord_flip() +
  scale_y_reverse() +
  geom_flow(aes(fill=group),width=.3, aes.flow = "backward") +
  geom_stratum(aes(fill=group),width=.3) +
  geom_text(stat = "stratum", aes(label = glue("{group}")),size=12) +
  ggrepel::geom_text_repel(data=plot_df_long_left,stat = "stratum", aes(label = glue("n={nn}\n({freq}%)")),nudge_x = -.3,
                           size=10,segment.color = NA) +
  ggrepel::geom_text_repel(data=plot_df_long_right,stat = "stratum", aes(label = glue("n={nn}\n({freq}%)")),nudge_x = .3,
                           size=10,segment.color = NA) +
  scale_fill_manual(name = "group",values=c("gray90","gray8 0","gray70",cluster_colors))+
  theme_void() +
  rremove("legend")

p1

tiff(glue("figures/cell_level_alluvial_plot_downsampled_{n_cells}.tiff"), width=500,height=300, units = "mm", res=600)
print(p1)
dev.off()
# 
# alluvial_plot_df <- plot_df
# ################################################################################
# 
# ################################################################################
# 
# treatment_table <- alluvial_plot_df %>% 
#   filter(treatment_status %in% c("Naive","CTX ± ICI")) %>% 
#   count(subtype,treatment_status) %>% 
#   arrange(subtype) %>% 
#   pivot_wider(names_from = treatment_status,values_from = n) %>% 
#   column_to_rownames("subtype") %>% 
#   select(CTX ± ICI,Naive)
# 
# subtypes <- rownames(treatment_table)
# 
# results_list <- list()
# i <- 1
# for(i in 1:4){
#   
#   contin_table <- rbind(treatment_table[i,],colSums(treatment_table[-i,]))
#   
#   fisher_res <- fisher.test((contin_table+.5))
#   
#   results_list[["subtype"]] <- append(results_list[["subtype"]], subtypes[i])
#   results_list[["or"]] <- append(results_list[["or"]], fisher_res$estimate)
#   results_list[["pval"]] <- append(results_list[["pval"]], fisher_res$p.value)
#   results_list[["up_or"]] <- append(results_list[["up_or"]], fisher_res$conf.int[1])
#   results_list[["low_or"]] <- append(results_list[["low_or"]], fisher_res$conf.int[2])
#   
# }
# 
# 
# plot_df <- as.data.frame(results_list)
# 
# plot_df <- plot_df %>% 
#   mutate(padj = p.adjust(pval), method = "BH") %>% 
#   mutate(signif = ifelse(padj < 0.05, "s","ns")) %>% 
#   mutate(log_or = log(or)) %>% 
#   mutate(log_upper_or = log(up_or)) %>% 
#   mutate(log_lower_or = log(low_or)) %>% 
#   mutate(star_height = ifelse(log_or > 0, log_or+.1,log_or-.1))
# 
# plot_df$subtype <- factor(plot_df$subtype, levels=c("A","N","P",'I'))
# 
# p2 <- ggplot(plot_df,aes(x=log_or,y=fct_rev(subtype),color=subtype))+
#   geom_point(aes(shape = factor(signif)),size=12,fill="white",show.legend = F, stroke=3)+
#   scale_shape_manual(values = c("ns" = 1, "s" = 16)) +
#   geom_errorbarh(aes(xmin = log_lower_or, xmax = log_upper_or), height = 0.1, linewidth = 1,show.legend = F)+
#   geom_vline(xintercept = 0, linetype = 2)+
#   scale_color_manual(values = cluster_colors)+
#   xlim(-2,2)+
#   labs(y="Subtype",
#        x="log(OR)")+
#   theme_classic()+
#   annotate("text", x=-.85, y=4.5, label = "Naive", angle=0,size=8) +
#   annotate("text", x=.85, y=4.5, label = "CTX ± ICI", angle=0,size=8) +
#   theme(axis.text = element_text(size=22,angle = 0, hjust = 1),
#         axis.title = element_text(size=24),
#         axis.text.x = element_text(angle = 0, hjust = .5))
# 
# p2
# 
# sprintf("%.4f", plot_df$padj)
# ################################################################################
# 
# tarla_table <- alluvial_plot_df %>% 
#   filter(!is.na(tarla)) %>% 
#   count(subtype,tarla) %>% 
#   arrange(subtype) %>% 
#   pivot_wider(names_from = tarla,values_from = n) %>% 
#   column_to_rownames("subtype") %>% 
#   select(post,pre)
# 
# subtypes <- rownames(tarla_table)
# 
# results_list <- list()
# i <- 1
# for(i in 1:4){
#   
#   contin_table <- rbind(tarla_table[i,],colSums(tarla_table[-i,]))
#   
#   fisher_res <- fisher.test((contin_table+.5))
# 
#   results_list[["subtype"]] <- append(results_list[["subtype"]], subtypes[i])
#   results_list[["or"]] <- append(results_list[["or"]], fisher_res$estimate)
#   results_list[["pval"]] <- append(results_list[["pval"]], fisher_res$p.value)
#   results_list[["up_or"]] <- append(results_list[["up_or"]], fisher_res$conf.int[1])
#   results_list[["low_or"]] <- append(results_list[["low_or"]], fisher_res$conf.int[2])
#   
# }
# 
# 
# plot_df <- as.data.frame(results_list)
# 
# plot_df <- plot_df %>% 
#   mutate(padj = p.adjust(pval), method = "BH") %>% 
#   mutate(signif = ifelse(padj < 0.05, "s","ns")) %>% 
#   mutate(log_or = log(or)) %>% 
#   mutate(log_upper_or = log(up_or)) %>% 
#   mutate(log_lower_or = log(low_or)) %>% 
#   mutate(star_height = ifelse(log_or > 0, log_or+.1,log_or-.1))
# 
# plot_df$subtype <- factor(plot_df$subtype, levels=c("A","N","P",'I'))
# 
# p3 <- ggplot(plot_df,aes(x=log_or,y=fct_rev(subtype),color=subtype))+
#   geom_point(aes(shape = factor(signif)),size=12,fill="white",show.legend = F, stroke=3)+
#   scale_shape_manual(values = c("ns" = 1, "s" = 16)) +
#   geom_errorbarh(aes(xmin = log_lower_or, xmax = log_upper_or), height = 0.1, linewidth = 1,show.legend = F)+
#   geom_vline(xintercept = 0, linetype = 2)+
#   scale_color_manual(values = cluster_colors)+
#   xlim(-2,2)+
#   labs(y="Subtype",
#        x="log(OR)")+
#   theme_classic()+
#   annotate("text", x=-1, y=4.35, label = "Pre\ntarlatamab", angle=0,size=8) +
#   annotate("text", x=1, y=4.35, label = "Post\ntarlatamab", angle=0,size=8) +
#   theme(axis.text = element_text(size=22,angle = 0, hjust = 1),
#         axis.title = element_text(size=24),
#         axis.text.x = element_text(angle = 0, hjust = .5))
# 
# p3
# 
# sprintf("%.4f", plot_df$padj)
# 
# tiff(glue("figures/downsampled_subtype_treatment_status_fisher_or_results_{n_cells}.tiff"), width=120,height=200, units = "mm", res=600)
# print(p2)
# dev.off()
# 
# 
# tiff(glue("figures/downsampled_subtype_tarla_status_fisher_or_results_{n_cells}.tiff"), width=120,height=200, units = "mm", res=600)
# print(p3)
# dev.off()
# 
# ################################################################################
# # Count patients
# sampled_data$patient_id %>% 
#   unique() %>% 
#   length()
# 
# # Count LBs
# sampled_data$collection_id %>% 
#   unique() %>% 
#   length()
# 
# # Count cells
# nrow(sampled_data)
# 
# # Naive patients
# sampled_data %>% 
#   filter(treatment_status == "Naive") %>% 
#   pull(patient_id) %>% 
#   unique() %>% 
#   length()
# 
# # Naive LBs
# sampled_data %>% 
#   filter(treatment_status == "Naive") %>% 
#   pull(collection_id) %>% 
#   unique() %>% 
#   length()
# 
# naive_lbs <- sampled_data %>% 
#   filter(treatment_status == "Naive") %>% 
#   pull(collection_id) %>% 
#   unique() 
# 
# # Naive cells
# sampled_data %>% 
#   filter(treatment_status == "Naive") %>% 
#   nrow()
# 
# # CTX ± ICI patients
# sampled_data %>% 
#   filter(treatment_status == "CTX ± ICI") %>% 
#   pull(patient_id) %>% 
#   unique() %>% 
#   length()
# 
# # CTX ± ICI LBs
# sampled_data %>% 
#   filter(treatment_status == "CTX ± ICI") %>% 
#   pull(collection_id) %>% 
#   unique() %>% 
#   length()
# 
# # CTX ± ICI cells
# sampled_data %>% 
#   filter(treatment_status == "CTX ± ICI") %>% 
#   nrow()
# 
# # Pre-tarla patients
# sampled_data %>% 
#   filter(tarla == "pre") %>% 
#   pull(patient_id) %>% 
#   unique() %>% 
#   length()
# 
# # Pre-tarla LBs
# sampled_data %>% 
#   filter(tarla == "pre") %>% 
#   pull(collection_id) %>% 
#   unique() %>% 
#   length()
# 
# # Pre-tarla cells
# sampled_data %>% 
#   filter(tarla == "pre") %>% 
#   nrow()
# 
# # Post-tarla patients
# sampled_data %>% 
#   filter(treatment_status == "Tarla") %>% 
#   pull(patient_id) %>% 
#   unique() %>% 
#   length()
# 
# # Post-tarla LBs
# sampled_data %>% 
#   filter(treatment_status == "Tarla") %>% 
#   pull(collection_id) %>% 
#   unique() %>% 
#   length()
# 
# # Post-tarla cells
# sampled_data %>% 
#   filter(treatment_status == "Tarla") %>% 
#   nrow()
# 
# 
# sampled_data %>% 
#   count(treatment_status,subtype) %>% 
#   group_by(treatment_status) %>% 
#   mutate(total = sum(n)) %>% 
#   mutate(freq = (n/total)*100) %>% 
#   select(treatment_status,subtype,freq)
# 
# 
# sampled_data %>% 
#   filter(collection_id %in% naive_lbs) %>% 
#   count(collection_id)
