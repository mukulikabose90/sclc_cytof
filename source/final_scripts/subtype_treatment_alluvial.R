source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

cluster_colors <- c("#dd4b33", "#D1DACF", "#A8DADC", "#457B9D")

################################################################################
# Sample level
################################################################################
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
  select(collection_id,treatment_status,tarla) %>%
  distinct() %>% 
  merge(.,sample_level_subtype,by="collection_id") %>% 
  select(subtype,treatment_status.x,tarla) %>% 
  rename("treatment_status" = "treatment_status.x") %>% 
  cbind(1) %>% 
  rename("n" = "1") 
  

plot_df$treatment_status <- ifelse(plot_df$treatment_status == "naive","Naive","SOC")
plot_df$tarla <- ifelse(plot_df$tarla == "pre","Pre-Tarla","Tarla")

plot_df$treatment_status <- ifelse(is.na(plot_df$tarla), plot_df$treatment_status, ifelse(plot_df$tarla == "Pre-Tarla","SOC","Tarla"))

plot_df$treatment_status <- factor(plot_df$treatment_status, levels=c("Naive","SOC","Pre-Tarla","Tarla"))

plot_df$subtype <- factor(plot_df$subtype, levels=c("A","N","P",'I'))



plot_df_long <- to_lodes_form(data.frame(plot_df),
                              key = "category", value = "group", id = "cohort",
                              axes = 1:2)

ggplot(data = plot_df_long,
       aes(x = category, stratum = group, alluvium = cohort, y = n)) +
  geom_flow(aes(fill=group),width=.25) +
  geom_stratum(aes(fill=group),width=.25) +
  geom_text(stat = "stratum", aes(label = group),size=6) +
  scale_fill_manual(name = "group",values=c(cluster_colors,"gray90","gray80","gray70"))+
  theme_void() +
  rremove("legend")

################################################################################
# Cell level
################################################################################

plot_df <- ctcs@colData %>% 
  as.data.frame() %>% 
  select(subtype,treatment_status,tarla) %>% 
  cbind(1) %>% 
  rename("n" = "1")

plot_df$treatment_status <- ifelse(plot_df$treatment_status == "naive","Naive","SOC")
plot_df$tarla <- ifelse(plot_df$tarla == "pre","Pre-Tarla","Tarla")

plot_df$treatment_status <- ifelse(is.na(plot_df$tarla), plot_df$treatment_status, ifelse(plot_df$tarla == "Pre-Tarla","SOC","Tarla"))

plot_df$treatment_status <- factor(plot_df$treatment_status, levels=c("Naive","SOC","Pre-Tarla","Tarla"))

plot_df$subtype <- factor(plot_df$subtype, levels=c("A","N","P",'I'))


plot_df_long <- to_lodes_form(data.frame(plot_df),
                              key = "category", value = "group", id = "cohort",
                              axes = 1:2)

ggplot(data = plot_df_long,
       aes(x = category, stratum = group, alluvium = cohort, y = n)) +
  geom_flow(aes(fill=group),width=.25) +
  geom_stratum(aes(fill=group),width=.25) +
  geom_text(stat = "stratum", aes(label = group),size=6) +
  scale_fill_manual(name = "group",values=c(cluster_colors,"gray90","gray80","gray70"))+
  theme_void() +
  rremove("legend")

