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

all_data <- ctcs@colData %>% 
  as.data.frame()

all_data$treatment_status <- ifelse(all_data$treatment_status == "naive","Naive","CTX ± ICI")

all_data$treatment_status <- ifelse(is.na(all_data$tarla), all_data$treatment_status,
                                        ifelse(all_data$tarla == "pre", all_data$treatment_status, "Tarla"))


plot_results <- list()
################################################################################
# Naive vs CTX + ICI Association
################################################################################
data_df <- all_data %>% 
  filter(treatment_status %in% c("Naive","CTX ± ICI")) %>% 
  select(subtype,treatment_status,patient_id)

# count cells and patients
data_df %>% 
  count(treatment_status)

data_df %>% 
  count(patient_id,treatment_status) %>% 
  count(treatment_status)

data_df$treatment_status <- factor(data_df$treatment_status, levels = c("Naive","CTX ± ICI"))

results_list <- list()
for(curr_subtype in c("A","N","P","Mes")){
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
  
  res <- data.frame("subtype"=curr_subtype,"or"=or,"pval"=curr_pval,"up_or"=upper_or,"low_or"=lower_or)
  
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

plot_df$subtype <- factor(plot_df$subtype,levels=c("A","N","P","Mes"))

plot_df$comparison <- "Naive_CTX ± ICI"

plot_results <- append(plot_results,list(plot_df))


sprintf("%.4f", plot_df$padj)

################################################################################
# Naive vs Tarla  Association
################################################################################
data_df <- all_data %>% 
  filter(treatment_status %in% c("Naive","Tarla")) %>% 
  select(subtype,treatment_status,patient_id)

# count cells and patients
data_df %>% 
  count(treatment_status)

data_df %>% 
  count(patient_id,treatment_status) %>% 
  count(treatment_status)

data_df$treatment_status <- factor(data_df$treatment_status, levels = c("Naive","Tarla"))

results_list <- list()
for(curr_subtype in c("A","N","P","Mes")){
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

plot_df$subtype <- factor(plot_df$subtype,levels=c("A","N","P","Mes"))

plot_df$comparison <- "Naive_Tarla"

plot_results <- append(plot_results,list(plot_df))

sprintf("%.4f", plot_df$padj)
################################################################################
# Naive vs Tarla  Association
################################################################################
data_df <- all_data %>% 
  filter(treatment_status %in% c("CTX ± ICI","Tarla")) %>% 
  select(subtype,treatment_status,patient_id)

# count cells and patients
data_df %>% 
  count(treatment_status)

data_df %>% 
  count(patient_id,treatment_status) %>% 
  count(treatment_status)

data_df$treatment_status <- factor(data_df$treatment_status, levels = c("CTX ± ICI","Tarla"))

results_list <- list()
for(curr_subtype in c("A","N","P","Mes")){
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

plot_df$subtype <- factor(plot_df$subtype,levels=c("A","N","P","Mes"))


plot_df$comparison <- "CTX ± ICI_Tarla"

plot_results <- append(plot_results,list(plot_df))


sprintf("%.4f", plot_df$padj)


plot_df <- as.data.frame(do.call(rbind,plot_results))

plot_df <- plot_df %>% 
  mutate(padj = p.adjust(pval), method = "BH") %>% 
  mutate(signif = ifelse(padj < 0.05, "s","ns"))

plot_df$left_label <- as.character(sapply(plot_df$comparison, function(x) strsplit(x,"_")[[1]][1]))
plot_df$right_label <- as.character(sapply(plot_df$comparison, function(x) strsplit(x,"_")[[1]][2]))

plot_df$subtype <- factor(plot_df$subtype, levels=c("A","N","P",'Mes'))
plot_df$comparison <- factor(plot_df$comparison, levels = c("Naive_CTX ± ICI","Naive_Tarla","CTX ± ICI_Tarla"))

p1 <- ggplot(plot_df,aes(x=log_or,y=fct_rev(subtype),color=subtype))+
  geom_point(aes(shape = factor(signif)),size=12,fill="white",show.legend = F, stroke=3)+
  facet_wrap(~comparison)+
  scale_shape_manual(values = c("ns" = 1, "s" = 16)) +
  geom_errorbarh(aes(xmin = log_lower_or, xmax = log_upper_or), height = 0.1, linewidth = 1,show.legend = F)+
  geom_vline(xintercept = 0, linetype = 2)+
  scale_color_manual(values = cluster_colors)+
  xlim(-3,3)+
  labs(y="Subtype",
       x="log(OR)")+
  guides(color="none")+
  theme_classic()+
  geom_text(x=-1.5, y=4.35, aes(label = left_label), angle=0,size=8, color="black")+
  geom_text(x=1.5, y=4.35, aes(label = right_label), angle=0,size=8, color="black") +
  theme(axis.text = element_text(size=22,angle = 0, hjust = 1),
        axis.title = element_text(size=24),
        axis.text.x = element_text(angle = 0, hjust = .5),
        strip.background = element_blank(),
        strip.text = element_blank())


p1
tiff(glue("figures/subtype_treatment_status_lme_results_all_comp.tiff"), width=360,height=200, units = "mm", res=600)
print(p1)
dev.off()



