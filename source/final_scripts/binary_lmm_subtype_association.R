source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)

cluster_colors <- c("#dd4b33", "#F1FAEE", "#A8DADC", "#457B9D")
################################################################################
ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

################################################################################

sce <- ctcs


ctcs@colData %>% 
  as.data.frame() %>% 
  count(patient_id)

# Sex
data_df <- sce@colData %>% 
  as.data.frame() %>% 
  select(subtype,sex,patient_id) %>% 
  filter(!is.na(sex))


# Age
data_df <- sce@colData %>% 
  as.data.frame() %>% 
  select(subtype,age,patient_id) %>% 
  filter(!is.na(age)) %>% 
  mutate(age = as.numeric(as.character(age)))

# Treatment Status
data_df <- sce@colData %>% 
  as.data.frame() %>% 
  select(subtype,treatment_status,patient_id) %>% 
  filter(!is.na(treatment_status))

data_df$treatment_status <- factor(data_df$treatment_status, levels = c("naive","treated"))

# Tarla Status
data_df <- sce@colData %>% 
  as.data.frame() %>% 
  select(subtype,tarla,patient_id) %>% 
  filter(!is.na(tarla))

data_df$tarla <- factor(data_df$tarla, levels = c("pre","post"))



results_list <- list()
for(curr_subtype in c("A","N","P","I")){
  data_df$curr_subtype <- as.factor(as.integer(data_df$subtype == curr_subtype))
  
  formula_str <- glue("curr_subtype ~ {colnames(data_df)[2]} + (1 | patient_id)")
  # formula_str <- glue("curr_subtype ~ {group}")
  
  model <- glmer(
    formula = as.formula(formula_str),
    family = binomial(link = "logit"),
    data = data_df)
  
  # model <- glm(
  #   formula = as.formula(formula_str),
  #   family = binomial(link = "logit"),
  #   data = data_df)
  
  summary(model)
  
  or <- exp(fixef(model)[2])
  
  tidy_out <- tidy(model,effects='fixed')
  curr_pval <- tidy_out$p.value[2]
  
  
  lower_or <- exp(confint(model, parm = "beta_", method = "Wald")[2,][1])
  upper_or <- exp(confint(model, parm = "beta_", method = "Wald")[2,][2])
  
  res <- data.frame("subtype"=curr_subtype,"or"=or,"pval"=curr_pval,up_or=upper_or,low_or=lower_or)
  
  results_list <- append(results_list, list(res))
  
}

# Combine into one data frame
all_results <- bind_rows(results_list)


plot_df <- all_results %>% 
  mutate(signif = ifelse(pval < 0.05, 16,1)) %>% 
  mutate(log_or = log(or)) %>% 
  mutate(log_upper_or = log(up_or)) %>% 
  mutate(log_lower_or = log(low_or)) %>% 
  mutate(star_height = ifelse(log_or > 0, log_or+.1,log_or-.1))



ggplot(plot_df,aes(x=log_or,y=subtype))+
  geom_errorbarh(aes(xmin = log_lower_or, xmax = log_upper_or), height = 0.1, linewidth = .5)+
  geom_point(aes(x=log_or,y=subtype, shape = factor(signif),color=subtype),size=4,fill="white",show.legend = F)+
  scale_shape_manual(values = c("1" = 1, "16" = 16)) +
  geom_vline(xintercept = 0, linetype = 2)+
  xlim(-6,6)+
  scale_y_discrete(limits=rev)+
  # scale_color_manual(values = cluster_colors)+
  labs(y="Cluster",
       x="log(OR)")+
  theme_classic()+
  theme(axis.text = element_text(size=12,angle = 0, hjust = 1),
        axis.title = element_text(size=14),
        axis.text.x = element_text(angle = 0, hjust = .5))+
  ggtitle("Post/Pre")


p1







plot_df <- do.call(rbind,final_df) %>% 
  as.data.frame()

colnames(plot_df) <- c("subtype","group","prob","pval")
plot_df$prob <- as.numeric(plot_df$prob)
plot_df$pval <- as.numeric(plot_df$pval)


plot_df$subtype <- factor(plot_df$subtype, levels=rev(c("A","N","P","I")))

ggplot(plot_df)+
  geom_point(aes(x=group,y=subtype,size=prob,fill=pval), shape = 21, color = "black", stroke = 0.6) +
  scale_fill_gradient2(low = "#457B9D", mid = "white", high = "#dd4b33", midpoint = 0) +
  scale_size(range = c(2, 8)) +
  theme_minimal() +
  theme(axis.text = element_text(size=12,angle = 0, hjust = 1),
        axis.title = element_text(size=14),
        axis.text.x = element_text(size=10,angle = 0, hjust = .5))+
  labs(title = "",
       size = "Probability", fill = "P value",
       y = "Subtype",
       x = "")



