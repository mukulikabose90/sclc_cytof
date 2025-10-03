################################################################################
# This script tests the subtype association with multiple variables and plots 
# the results as log odds ratio forest plots
################################################################################
source("source/sclc_cytof_functions.R")

set.seed(42)

cluster_colors <- c("#dd4b33", "#D1DACF", "#A8DADC", "#457B9D")
################################################################################
#Read in data
################################################################################
ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

################################################################################
# Sex Association
################################################################################
data_df <- ctcs@colData %>% 
  as.data.frame() %>% 
  select(subtype,sex,patient_id) %>% 
  filter(!is.na(sex))

results_list <- list()
for(curr_subtype in c("A","N","P","I")){
  data_df$curr_subtype <- as.factor(as.integer(data_df$subtype == curr_subtype))
  
  formula_str <- glue("curr_subtype ~ {colnames(data_df)[2]} + (1 | patient_id)")
  
  model <- glmer(
    formula = as.formula(formula_str),
    family = binomial(link = "logit"),
    data = data_df)
  
  or <- exp(fixef(model)[2])
  
  tidy_out <- broom.mixed::tidy(model,effects='fixed')
  curr_pval <- tidy_out$p.value[2]
  
  lower_or <- exp(confint(model, parm = "beta_", method = "Wald")[2,][1])
  upper_or <- exp(confint(model, parm = "beta_", method = "Wald")[2,][2])
  
  res <- data.frame("subtype"=curr_subtype,"or"=or,"pval"=curr_pval,up_or=upper_or,low_or=lower_or)
  
  results_list <- append(results_list, list(res))
  
}

# Combine into one data frame
all_results <- bind_rows(results_list)

plot_df <- all_results %>% 
  mutate(padj = p.adjust(pval), method = "BH") %>% 
  mutate(signif = ifelse(padj < 0.05, "s","ns")) %>% 
  mutate(log_or = log(or)) %>% 
  mutate(log_upper_or = log(up_or)) %>% 
  mutate(log_lower_or = log(low_or)) %>% 
  mutate(star_height = ifelse(log_or > 0, log_or+.1,log_or-.1))

plot_df$subtype <- factor(plot_df$subtype,levels=c("A","N","P","I"))

p1 <- ggplot(plot_df,aes(x=log_or,y=fct_rev(subtype),color=subtype))+
  geom_point(aes(shape = factor(signif)),size=9,fill="white",show.legend = F)+
  scale_shape_manual(values = c("ns" = 1, "s" = 16)) +
  geom_errorbarh(aes(xmin = log_lower_or, xmax = log_upper_or), height = 0.1, linewidth = .5,show.legend = F)+
  geom_vline(xintercept = 0, linetype = 2)+
  scale_color_manual(values = cluster_colors)+
  xlim(-1.5,1.5)+
  labs(y="Subtype",
       x="log(OR)")+
  theme_classic()+
  annotate("text", x=-1, y=4.5, label = "Female", angle=0,size=5) +
  annotate("text", x=1, y=4.5, label = "Male", angle=0,size=5) +
  theme(axis.text = element_text(size=18,angle = 0, hjust = 1),
        axis.title = element_text(size=20),
        axis.text.x = element_text(angle = 0, hjust = .5))

tiff("figures/subtype_sex_or_results.tiff", width=100,height=200, units = "mm", res=600)
print(p1)
dev.off()
################################################################################
# Age Association
################################################################################
data_df <- ctcs@colData %>% 
  as.data.frame() %>% 
  select(subtype,age,patient_id) %>% 
  filter(!is.na(age)) %>% 
  mutate(age = as.numeric(as.character(age)))

results_list <- list()
for(curr_subtype in c("A","N","P","I")){
  data_df$curr_subtype <- as.factor(as.integer(data_df$subtype == curr_subtype))
  
  formula_str <- glue("curr_subtype ~ {colnames(data_df)[2]} + (1 | patient_id)")
  
  model <- glmer(
    formula = as.formula(formula_str),
    family = binomial(link = "logit"),
    data = data_df)
  
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
  mutate(padj = p.adjust(pval), method = "BH") %>% 
  mutate(signif = ifelse(padj < 0.05, "s","ns")) %>% 
  mutate(log_or = log(or)) %>% 
  mutate(log_upper_or = log(up_or)) %>% 
  mutate(log_lower_or = log(low_or)) %>% 
  mutate(star_height = ifelse(log_or > 0, log_or+.1,log_or-.1))

plot_df$subtype <- factor(plot_df$subtype,levels=c("A","N","P","I"))

p2 <- ggplot(plot_df,aes(x=log_or,y=fct_rev(subtype),color=subtype))+
  geom_point(aes(shape = factor(signif)),size=9,fill="white",show.legend = F)+
  scale_shape_manual(values = c("ns" = 1, "s" = 16)) +
  geom_errorbarh(aes(xmin = log_lower_or, xmax = log_upper_or), height = 0.1, linewidth = .5,show.legend = F)+
  geom_vline(xintercept = 0, linetype = 2)+
  scale_color_manual(values = cluster_colors)+
  xlim(-1,1)+
  labs(y="Subtype",
       x="log(OR)")+
  theme_classic()+
  annotate("text", x=-.5, y=4.5, label = "Younger", angle=0,size=5) +
  annotate("text", x=.5, y=4.5, label = "Older", angle=0,size=5) +
  theme(axis.text = element_text(size=18,angle = 0, hjust = 1),
        axis.title = element_text(size=20),
        axis.text.x = element_text(angle = 0, hjust = .5))

tiff("figures/subtype_age_or_results.tiff", width=100,height=200, units = "mm", res=600)
print(p2)
dev.off()

################################################################################
# Treatment Status Association
################################################################################
data_df <- ctcs@colData %>% 
  as.data.frame() %>% 
  select(subtype,treatment_status,patient_id,tarla) %>% 
  filter(tarla == "pre" | is.na(tarla))

# count cells and patients
data_df %>% 
  count(treatment_status)

data_df %>% 
  count(patient_id,treatment_status) %>% 
  count(treatment_status)

data_df$treatment_status <- factor(data_df$treatment_status, levels = c("naive","treated"))

results_list <- list()
for(curr_subtype in c("A","N","P","I")){
  data_df$curr_subtype <- as.factor(as.integer(data_df$subtype == curr_subtype))
  
  formula_str <- glue("curr_subtype ~ {colnames(data_df)[2]} + (1 | patient_id)")

  model <- glmer(
    formula = as.formula(formula_str),
    family = binomial(link = "logit"),
    data = data_df)
  
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
  mutate(padj = p.adjust(pval), method = "BH") %>% 
  mutate(signif = ifelse(padj < 0.05, "s","ns")) %>% 
  mutate(log_or = log(or)) %>% 
  mutate(log_upper_or = log(up_or)) %>% 
  mutate(log_lower_or = log(low_or)) %>% 
  mutate(star_height = ifelse(log_or > 0, log_or+.1,log_or-.1))

plot_df$subtype <- factor(plot_df$subtype,levels=c("A","N","P","I"))

p3 <- ggplot(plot_df,aes(x=log_or,y=fct_rev(subtype),color=subtype))+
  geom_point(aes(shape = factor(signif)),size=12,fill="white",show.legend = F, stroke=3)+
  scale_shape_manual(values = c("ns" = 1, "s" = 16)) +
  geom_errorbarh(aes(xmin = log_lower_or, xmax = log_upper_or), height = 0.1, linewidth = 1,show.legend = F)+
  geom_vline(xintercept = 0, linetype = 2)+
  scale_color_manual(values = cluster_colors)+
  xlim(-2,2)+
  labs(y="Subtype",
       x="log(OR)")+
  theme_classic()+
  annotate("text", x=-.85, y=4.5, label = "Naive", angle=0,size=6) +
  annotate("text", x=.85, y=4.5, label = "SOC", angle=0,size=6) +
  theme(axis.text = element_text(size=22,angle = 0, hjust = 1),
        axis.title = element_text(size=24),
        axis.text.x = element_text(angle = 0, hjust = .5))


tiff("figures/subtype_treatment_status_or_results.tiff", width=100,height=200, units = "mm", res=600)
print(p3)
dev.off()

################################################################################
# Tarla Status Association
################################################################################
data_df <- ctcs@colData %>% 
  as.data.frame() %>% 
  select(subtype,tarla,patient_id) %>% 
  filter(!is.na(tarla))

data_df$tarla <- factor(data_df$tarla, levels = c("pre","post"))

# count cells and patients
data_df %>% 
  count(tarla)

data_df %>% 
  count(patient_id,tarla) %>% 
  count(tarla)

results_list <- list()
for(curr_subtype in c("A","N","P","I")){
  data_df$curr_subtype <- as.factor(as.integer(data_df$subtype == curr_subtype))
  
  formula_str <- glue("curr_subtype ~ {colnames(data_df)[2]} + (1 | patient_id)")
  
  model <- glmer(
    formula = as.formula(formula_str),
    family = binomial(link = "logit"),
    data = data_df)
  
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
  mutate(padj = p.adjust(pval), method = "BH") %>% 
  mutate(signif = ifelse(padj < 0.05, "s","ns")) %>% 
  mutate(log_or = log(or)) %>% 
  mutate(log_upper_or = log(up_or)) %>% 
  mutate(log_lower_or = log(low_or)) %>% 
  mutate(star_height = ifelse(log_or > 0, log_or+.1,log_or-.1))

plot_df$subtype <- factor(plot_df$subtype,levels=c("A","N","P","I"))

p4 <- ggplot(plot_df,aes(x=log_or,y=fct_rev(subtype),color=subtype))+
  geom_point(aes(shape = factor(signif)),size=12,fill="white",show.legend = F,stroke=3)+
  scale_shape_manual(values = c("ns" = 1, "s" = 16)) +
  geom_errorbarh(aes(xmin = log_lower_or, xmax = log_upper_or), height = 0.1, linewidth = 1,show.legend = F)+
  geom_vline(xintercept = 0, linetype = 2)+
  scale_color_manual(values = cluster_colors)+
  xlim(-1.5,1.5)+
  labs(y="Subtype",
       x="log(OR)")+
  theme_classic()+
  annotate("text", x=-.75, y=4.5, label = "Pre-Tarlatamab", angle=0,size=4) +
  annotate("text", x=.75, y=4.5, label = "Post-Tarlatamab", angle=0,size=4) +
  theme(axis.text = element_text(size=22,angle = 0, hjust = 1),
        axis.title = element_text(size=24),
        axis.text.x = element_text(angle = 0, hjust = .5)) 
  
  
tiff("figures/subtype_tarla_or_results.tiff", width=100, height=200, units = "mm", res=600)
print(p4)
dev.off()



