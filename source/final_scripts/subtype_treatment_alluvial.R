################################################################################
# This script plots alluvial plots indicating the proportions of subtypes and treatment status.
#   1) Cell level alluvial
#   2) Sample level alluvial
#   3) Sample level alluvial (excluding longitudinal patients)
################################################################################
source("source/sclc_cytof_functions.R")

set.seed(42)
################################################################################
# Read in data
################################################################################
ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

cluster_colors <- c("#dd4b33", "#D1DACF", "#A8DADC", "#457B9D")
################################################################################
# Cell level alluvial
################################################################################

ctcs <- ctcs[,ctcs$collection_id != "MDA-SC399-2"]


plot_df <- ctcs@colData %>% 
  as.data.frame() %>% 
  select(treatment_status,subtype,tarla) %>% 
  cbind(1) %>% 
  rename("n" = "1")

plot_df$treatment_status <- ifelse(plot_df$treatment_status == "naive","Naive","CTX ± ICI")

plot_df$treatment_status <- ifelse(is.na(plot_df$tarla), plot_df$treatment_status,
                                        ifelse(plot_df$tarla == "pre", plot_df$treatment_status, "Tarla"))

plot_df$subtype <- factor(plot_df$subtype, levels=c("A","N","P","Mes"))

# Count cells
plot_df %>% 
  count(treatment_status)

plot_df %>% 
  count(treatment_status,subtype) %>% 
  group_by(treatment_status) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n/total)*100) %>% 
  select(treatment_status,subtype,freq)

plot_df_long <- to_lodes_form(data.frame(plot_df),
                              key = "category", value = "group", id = "cohort",
                              axes = 1:2) %>% 
  add_count(group) %>% 
  group_by(category) %>% 
  mutate(total=sum(n)) %>% 
  mutate(freq=sprintf("%.1f",(nn/total)*100)) 



plot_df_long_left <- subset(plot_df_long, group %in% c("Naive","CTX ± ICI","Tarla"))
plot_df_long_right <- subset(plot_df_long, !group %in% c("Naive","CTX ± ICI","Tarla"))


plot_df_long$group <- factor(plot_df_long$group, levels=c("Naive","CTX ± ICI","Tarla", "A", "N", "P", "Mes"))
plot_df_long_left$group <- factor(plot_df_long_left$group, levels=c("Naive","CTX ± ICI","Tarla", "A", "N", "P", "Mes"))
plot_df_long_right$group <- factor(plot_df_long_right$group, levels=c("Naive","CTX ± ICI","Tarla", "A", "N", "P", "Mes"))

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

tiff("figures/cell_level_alluvial_plot.tiff", width=400,height=250, units = "mm", res=600)
print(p1)
dev.off()
################################################################################
# Sample level alluvial
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
  select(collection_id,treatment_status,tarla,patient_id) %>%
  distinct() %>% 
  merge(.,sample_level_subtype,by="collection_id") %>% 
  select(treatment_status.x,subtype,tarla,patient_id) %>% 
  rename("treatment_status" = "treatment_status.x") %>% 
  cbind(1) %>% 
  rename("n" = "1") 


plot_df$treatment_status <- ifelse(plot_df$treatment_status == "naive","Naive","CTX ± ICI")

plot_df$treatment_status <- ifelse(is.na(plot_df$tarla), plot_df$treatment_status,
                                   ifelse(plot_df$tarla == "pre", plot_df$treatment_status, "Tarla"))

plot_df$subtype <- factor(plot_df$subtype, levels=c("A","N","P","Mes"))

# Count cells
plot_df %>% 
  count(treatment_status)

plot_df_long <- to_lodes_form(data.frame(plot_df),
                              key = "category", value = "group", id = "cohort",
                              axes = 1:2) %>% 
  add_count(group) %>% 
  group_by(category) %>% 
  mutate(total=sum(n)) %>% 
  mutate(freq=sprintf("%.1f",(nn/total)*100)) %>% 
  mutate(stats = ifelse(group == "A", "", glue("\U00A0\n\n\nn={nn} ({freq}%)"))) %>% 
  mutate(a_stats = ifelse(group == "A",glue("n={nn} ({freq}%)"),""))



plot_df_long_left <- subset(plot_df_long, group %in% c("Naive","CTX ± ICI","Tarla"))
plot_df_long_right <- subset(plot_df_long, !group %in% c("Naive","CTX ± ICI","Tarla"))

plot_df_long$group <- factor(plot_df_long$group, levels=c("Naive","CTX ± ICI","Tarla", "A", "N", "P", "Mes"))
plot_df_long_left$group <- factor(plot_df_long_left$group, levels=c("Naive","CTX ± ICI","Tarla", "A", "N", "P", "Mes"))
plot_df_long_right$group <- factor(plot_df_long_right$group, levels=c("Naive","CTX ± ICI","Tarla", "A", "N", "P", "Mes"))

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
# Sample level (non-longitudinal patients) alluvial
################################################################################
# Find patient level subtypes
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

# Count cells
samples_used <- as.character(patient_level_subtype$collection_id)

count_df <- ctcs@colData %>% 
  as.data.frame() %>% 
  filter(collection_id %in% samples_used) 

count_df$treatment_status <- ifelse(count_df$treatment_status == "naive","Naive","CTX ± ICI")
count_df$tarla <- ifelse(count_df$tarla == "pre","Pre-Tarla","Tarla")

count_df$treatment_status <- ifelse(count_df$treatment_status == "Naive", count_df$treatment_status, ifelse(count_df$tarla == "Pre-Tarla" | is.na(count_df$tarla),count_df$treatment_status,"Tarla"))

count_df$treatment_status <- factor(count_df$treatment_status, levels=c("Naive","CTX ± ICI","Pre-Tarla","Tarla"))

count_df %>% 
  count(treatment_status)

count_df %>% 
  count(subtype)

#######

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


plot_df$treatment_status <- ifelse(plot_df$treatment_status == "naive","Naive","CTX ± ICI")

plot_df$treatment_status <- ifelse(is.na(plot_df$tarla), plot_df$treatment_status,
                                   ifelse(plot_df$tarla == "pre", plot_df$treatment_status, "Tarla"))

plot_df$subtype <- factor(plot_df$subtype, levels=c("A","N","P","Mes"))


plot_df_long <- to_lodes_form(data.frame(plot_df),
                              key = "category", value = "group", id = "cohort",
                              axes = 1:2) %>% 
  add_count(group) %>% 
  group_by(category) %>% 
  mutate(total=sum(n)) %>% 
  mutate(freq=sprintf("%.1f",(nn/total)*100)) 



plot_df_long_left <- subset(plot_df_long, group %in% c("Naive","CTX ± ICI","Tarla"))
plot_df_long_right <- subset(plot_df_long, !group %in% c("Naive","CTX ± ICI","Tarla"))


p3 <- ggplot(data = plot_df_long,
             aes(x = category, stratum = group, alluvium = cohort, y = total)) +
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
tiff("figures/patient_level_alluvial_plot.tiff", width=300,height=500, units = "mm", res=600)
print(p3)
dev.off()

