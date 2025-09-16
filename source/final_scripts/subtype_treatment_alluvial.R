source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

cluster_colors <- c("#dd4b33", "#D1DACF", "#A8DADC", "#457B9D")
################################################################################
# Cell level
################################################################################
ctcs@colData %>% 
  as.data.frame() %>% 
  pull(patient_id) %>% 
  unique() %>% 
  length()



plot_df <- ctcs@colData %>% 
  as.data.frame() %>% 
  select(treatment_status,subtype,tarla) %>% 
  cbind(1) %>% 
  rename("n" = "1")

plot_df$treatment_status <- ifelse(plot_df$treatment_status == "naive","Naive","SOC")
plot_df$tarla <- ifelse(plot_df$tarla == "pre","Pre-Tarla","Tarla")

plot_df$treatment_status <- ifelse(plot_df$treatment_status == "Naive", plot_df$treatment_status, ifelse(plot_df$tarla == "Pre-Tarla" | is.na(plot_df$tarla),plot_df$treatment_status,"Tarla"))

plot_df$treatment_status <- factor(plot_df$treatment_status, levels=c("Naive","SOC","Pre-Tarla","Tarla"))

plot_df$subtype <- factor(plot_df$subtype, levels=c("A","N","P",'I'))


# to_lodes_form(data.frame(plot_df),
#               key = "category", value = "group", id = "cohort",
#               axes = 1:2) %>% 
#   add_count(group)


plot_df_long <- to_lodes_form(data.frame(plot_df),
                              key = "category", value = "group", id = "cohort",
                              axes = 1:2) %>% 
  add_count(group) %>% 
  group_by(category) %>% 
  mutate(total=sum(n)) %>% 
  mutate(freq=sprintf("%.1f",(nn/total)*100)) 



plot_df_long_left <- subset(plot_df_long, group %in% c("Naive","SOC","Tarla"))
plot_df_long_right <- subset(plot_df_long, !group %in% c("Naive","SOC","Tarla"))


# plot_df_long$group <- factor(plot_df_long$group, levels = c("Tarla","SOC","Naive","I","P","N","A"))

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

# "\U00A0\n\nn={nn}\n({freq}%)"

tiff("figures/cell_level_alluvial_plot.tiff", width=400,height=250, units = "mm", res=600)
print(p1)
dev.off()
################################################################################
# Sample level
################################################################################
samples_to_use <- ctcs@colData %>% 
  as.data.frame() %>% 
  count(collection_id) %>% 
  filter(n > 10) %>% 
  pull(collection_id) %>% 
  as.character() %>% 
  unique()


sample_level_subtype <- ctcs@colData %>% 
  as.data.frame() %>% 
  select(collection_id,treatment_status,subtype) %>% 
  group_by(treatment_status) %>% 
  count(collection_id,subtype) %>% 
  group_by(collection_id,treatment_status) %>% 
  filter(n == max(n)) %>% 
  select(collection_id,subtype) 



plot_df <- ctcs@colData %>% 
  as.data.frame() %>% 
  filter(collection_id %in% samples_to_use) %>% 
  select(collection_id,treatment_status,tarla) %>%
  distinct() %>% 
  merge(.,sample_level_subtype,by="collection_id") %>% 
  select(treatment_status.x,subtype,tarla) %>% 
  rename("treatment_status" = "treatment_status.x") %>% 
  cbind(1) %>% 
  rename("n" = "1") 


plot_df$treatment_status <- ifelse(plot_df$treatment_status == "naive","Naive","SOC")
plot_df$tarla <- ifelse(plot_df$tarla == "pre","Pre-Tarla","Tarla")

plot_df$treatment_status <- ifelse(plot_df$treatment_status == "Naive", plot_df$treatment_status, ifelse(plot_df$tarla == "Pre-Tarla" | is.na(plot_df$tarla),plot_df$treatment_status,"Tarla"))

plot_df$treatment_status <- factor(plot_df$treatment_status, levels=c("Naive","SOC","Pre-Tarla","Tarla"))

plot_df$subtype <- factor(plot_df$subtype, levels=c("A","N","P",'I'))



plot_df_long <- to_lodes_form(data.frame(plot_df),
                              key = "category", value = "group", id = "cohort",
                              axes = 1:2) %>% 
  add_count(group) %>% 
  group_by(category) %>% 
  mutate(total=sum(n)) %>% 
  mutate(freq=sprintf("%.1f",(nn/total)*100)) %>% 
  mutate(stats = ifelse(group == "A", "", glue("\U00A0\n\n\nn={nn} ({freq}%)"))) %>% 
  mutate(a_stats = ifelse(group == "A",glue("n={nn} ({freq}%)"),""))



plot_df_long_left <- subset(plot_df_long, group %in% c("Naive","SOC","Tarla"))
plot_df_long_right <- subset(plot_df_long, !group %in% c("Naive","SOC","Tarla"))


p2 <- ggplot(data = plot_df_long,
             aes(x = category, stratum = group, alluvium = cohort, y = total)) +
  coord_flip() +
  scale_y_reverse() +
  geom_flow(aes(fill=group),width=.3, aes.flow = "backward") +
  geom_stratum(aes(fill=group),width=.3) +
  geom_text(stat = "stratum", aes(label = glue("{group}")),size=15) +
  ggrepel::geom_text_repel(data=plot_df_long_left,stat = "stratum", aes(label = glue("n={nn}\n({freq}%)")),nudge_x = -.3,
                           size=10,segment.color = NA) +
  ggrepel::geom_text_repel(data=plot_df_long_right,stat = "stratum", aes(label = glue("n={nn}\n({freq}%)")),nudge_x = .3,
                           size=10,segment.color = NA) +
  scale_fill_manual(name = "group",values=c("gray90","gray8 0","gray70",cluster_colors))+
  theme_void() +
  rremove("legend")

p2

tiff("figures/sample_level_alluvial_plot.tiff", width=500,height=300, units = "mm", res=600)
print(p2)
dev.off()

################################################################################
# Sample level (non-longitudinal patient)
################################################################################
samples_to_use <- ctcs@colData %>% 
  as.data.frame() %>% 
  count(collection_id) %>% 
  filter(n > 10) %>% 
  pull(collection_id) %>% 
  as.character() %>% 
  unique()

non_long_patients <- ctcs@colData %>%
  as.data.frame() %>% 
  filter(collection_id %in% samples_to_use) %>% 
  select(patient_id,sample_num) %>% 
  distinct() %>% 
  count(patient_id) %>% 
  filter(n == 1) %>% 
  pull(patient_id) %>% 
  as.character() 


patient_level_subtype <- ctcs@colData %>% 
  as.data.frame() %>% 
  filter(patient_id %in% non_long_patients & collection_id %in% samples_to_use) %>% 
  select(collection_id,subtype) %>% 
  count(collection_id,subtype) %>% 
  group_by(collection_id) %>% 
  filter(n == max(n)) %>% 
  select(collection_id,subtype) %>% 
  distinct()


tied_patients <- patient_level_subtype %>% 
  count(collection_id) %>% 
  filter(n > 1)



plot_df <- ctcs@colData %>% 
  as.data.frame() %>% 
  select(collection_id,treatment_status,tarla) %>%
  distinct() %>% 
  merge(.,patient_level_subtype,by="collection_id") %>% 
  select(treatment_status,subtype,tarla) %>% 
  cbind(1) %>% 
  rename("n" = "1") 


plot_df$treatment_status <- ifelse(plot_df$treatment_status == "naive","Naive","SOC")
plot_df$tarla <- ifelse(plot_df$tarla == "pre","Pre-Tarla","Tarla")

plot_df$treatment_status <- ifelse(plot_df$treatment_status == "Naive", plot_df$treatment_status, ifelse(plot_df$tarla == "Pre-Tarla" | is.na(plot_df$tarla),plot_df$treatment_status,"Tarla"))

plot_df$treatment_status <- factor(plot_df$treatment_status, levels=c("Naive","SOC","Pre-Tarla","Tarla"))

plot_df$subtype <- factor(plot_df$subtype, levels=c("A","N","P",'I'))



plot_df_long <- to_lodes_form(data.frame(plot_df),
                              key = "category", value = "group", id = "cohort",
                              axes = 1:2) %>% 
  add_count(group) %>% 
  group_by(category) %>% 
  mutate(total=sum(n)) %>% 
  mutate(freq=sprintf("%.1f",(nn/total)*100)) 



plot_df_long_left <- subset(plot_df_long, group %in% c("Naive","SOC","Tarla"))
plot_df_long_right <- subset(plot_df_long, !group %in% c("Naive","SOC","Tarla"))


p3 <- ggplot(data = plot_df_long,
             aes(x = category, stratum = group, alluvium = cohort, y = total)) +
  coord_flip() +
  scale_y_reverse() +
  geom_flow(aes(fill=group),width=.3, aes.flow = "backward") +
  geom_stratum(aes(fill=group),width=.3) +
  geom_text(stat = "stratum", aes(label = glue("{group}")),size=15) +
  ggrepel::geom_text_repel(data=plot_df_long_left,stat = "stratum", aes(label = glue("n={nn}\n({freq}%)")),nudge_x = -.3,
                           size=10,segment.color = NA) +
  ggrepel::geom_text_repel(data=plot_df_long_right,stat = "stratum", aes(label = glue("n={nn}\n({freq}%)")),nudge_x = .3,
                           size=10,segment.color = NA) +
  scale_fill_manual(name = "group",values=c("gray90","gray8 0","gray70",cluster_colors))+
  theme_void() +
  rremove("legend")

p3

tiff("figures/patient_level_alluvial_plot.tiff", width=500,height=300, units = "mm", res=600)
print(p3)
dev.off()

