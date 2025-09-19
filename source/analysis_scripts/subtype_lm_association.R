source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)

# cluster_colors <- c("#dd4b33", "#F1FAEE", "#A8DADC", "#457B9D")
cluster_colors <- c("#dd4b33", "#D1DACF", "#A8DADC", "#457B9D")
################################################################################
ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")
################################################################################

sce <- ctcs

################################################################################
# Sex
################################################################################
data_df <- sce@colData %>% 
  as.data.frame() %>% 
  select(subtype,sex,sample_id) %>% 
  filter(!is.na(sex))

wide_data_df <- data_df %>% 
  count(subtype,sample_id) %>% 
  group_by(sample_id) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n/total)*100) %>% 
  select(subtype,sample_id,freq) %>% 
  pivot_wider(names_from = "subtype",values_from = "freq",values_fill = 0) %>% 
  merge(.,data_df,by='sample_id') %>% 
  select(sex,A,N,P,I,sample_id) %>% 
  distinct()

results_list <- list()
pvals <- c()
for(curr_subtype in c("A","N","P","I")){
  
  formula_str <- glue("{curr_subtype} ~ {colnames(wide_data_df)[1]}")
  
  model <- lm(formula_str, data = wide_data_df)
  
  model_summary <- summary(model)
  
  curr_pval <- model_summary$coefficients[2,4]
  
  
  pvals <- append(pvals,curr_pval)
  
  
  
  results_list <- append(results_list, list(model))
  
}


names(pvals) <- c("A","N","P","I")

pval_df <- pvals %>% 
  enframe() %>% 
  rename("subtype"="name","pval"="value")

plot_df <- wide_data_df %>% 
  pivot_longer(!c(sample_id,sex), names_to = "subtype",values_to = "freq") %>% 
  merge(.,pval_df,by="subtype")

plot_df$subtype <- factor(plot_df$subtype, levels = c("A","N","P","I"))

p1 <- ggplot(plot_df, aes(x = sex, y = freq)) +
  facet_wrap(~subtype,nrow=1)+
  geom_boxplot(aes(fill=sex),outlier.shape = NA) +
  geom_text(aes(label=glue("p-value = {sprintf('%.2f',pval)}"),x=1.5,y=103,size=5))+
  geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
  scale_fill_manual(values=c("#A8DADC", "#dd4b33"))+
  theme_classic() +
  labs(y = "Subtype %", 
       x = "Sex")+
  rremove("legend")+
  theme(axis.text = element_text(size=18,angle = 0, hjust = 1),
        axis.title = element_text(size=20),
        axis.text.x = element_text(angle = 0, hjust = .5),
        strip.text = element_text(size=20))

p1

tiff("figures/subtype_sex_percent_results.tiff", width=300,height=100, units = "mm", res=600)
print(p1)
dev.off()

################################################################################
# Age
################################################################################
data_df <- sce@colData %>% 
  as.data.frame() %>% 
  select(subtype,age,sample_id) %>% 
  filter(!is.na(age))

wide_data_df <- data_df %>% 
  count(subtype,sample_id) %>% 
  group_by(sample_id) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n/total)*100) %>% 
  select(subtype,sample_id,freq) %>% 
  pivot_wider(names_from = "subtype",values_from = "freq",values_fill = 0) %>% 
  merge(.,data_df,by='sample_id') %>% 
  select(age,A,N,P,I,sample_id) %>% 
  distinct()

 
wide_data_df$age <- as.numeric(as.character(wide_data_df$age))

summary(wide_data_df$age)

results_list <- list()
pvals <- c()
for(curr_subtype in c("A","N","P","I")){
  
  formula_str <- glue("{curr_subtype} ~ {colnames(wide_data_df)[1]}")
  
  model <- lm(formula_str, data = wide_data_df)
  
  model_summary <- summary(model)
  
  curr_pval <- model_summary$coefficients[2,4]
  
  
  pvals <- append(pvals,curr_pval)
  
  results_list <- append(results_list, list(model))
  
}

summary(wide_data_df$age)
names(pvals) <- c("A","N","P","I")

pval_df <- pvals %>% 
  enframe() %>% 
  rename("subtype"="name","pval"="value")

plot_df <- wide_data_df %>% 
  pivot_longer(!c(sample_id,age), names_to = "subtype",values_to = "freq") %>% 
  merge(.,pval_df,by="subtype")

plot_df$subtype <- factor(plot_df$subtype, levels = c("A","N","P","I"))

p2 <- ggplot(plot_df, aes(x = age, y = freq)) +
  facet_wrap(~subtype,nrow=1)+
  geom_point()+
  geom_smooth(method = "lm")+
  geom_text(aes(label=glue("p-value = {sprintf('%.2f',pval)}"),x=62.5,y=105,size=5))+
  geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
  scale_fill_manual(values=c("#A8DADC", "#dd4b33"))+
  theme_classic() +
  xlim(37,88)+
  labs(y = "Subtype %", 
       x = "Age")+
  rremove("legend")+
  theme(axis.text = element_text(size=18,angle = 0, hjust = 1),
        axis.title = element_text(size=20),
        axis.text.x = element_text(angle = 0, hjust = .5),
        strip.text = element_text(size=20))

p2

tiff("figures/subtype_age_percent_results.tiff", width=300,height=100, units = "mm", res=600)
print(p2)
dev.off()
 ################################################################################
# Treatment Status
################################################################################
data_df <- sce@colData %>% 
  as.data.frame() %>% 
  select(subtype,treatment_status,sample_id) %>% 
  filter(!is.na(treatment_status))

wide_data_df <- data_df %>% 
  count(subtype,sample_id) %>% 
  group_by(sample_id) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n/total)*100) %>% 
  select(subtype,sample_id,freq) %>% 
  pivot_wider(names_from = "subtype",values_from = "freq",values_fill = 0) %>% 
  merge(.,data_df,by='sample_id') %>% 
  select(treatment_status,A,N,P,I,sample_id) %>% 
  distinct()

results_list <- list()
pvals <- c()
for(curr_subtype in c("A","N","P","I")){
  
  formula_str <- glue("{curr_subtype} ~ {colnames(wide_data_df)[1]}")
  
  model <- lm(formula_str, data = wide_data_df)
  
  model_summary <- summary(model)
  
  curr_pval <- model_summary$coefficients[2,4]
  
  pvals <- append(pvals,curr_pval)
  
  results_list <- append(results_list, list(model))
  
}


names(pvals) <- c("A","N","P","I")

pval_df <- pvals %>% 
  enframe() %>% 
  rename("subtype"="name","pval"="value")

plot_df <- wide_data_df %>% 
  pivot_longer(!c(sample_id,treatment_status), names_to = "subtype",values_to = "freq") %>% 
  merge(.,pval_df,by="subtype")

plot_df$subtype <- factor(plot_df$subtype, levels = c("A","N","P","I"))
plot_df$treatment_status <- ifelse(plot_df$treatment_status == "naive","Naive","Treated")

p3 <- ggplot(plot_df, aes(x = treatment_status, y = freq)) +
  facet_wrap(~subtype,nrow=1)+
  geom_boxplot(aes(fill=treatment_status),outlier.shape = NA) +
  geom_text(aes(label=glue("p-value = {sprintf('%.2f',pval)}"),x=1.5,y=103,size=5))+
  geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
  scale_fill_manual(values=c("#A8DADC", "#dd4b33"))+
  theme_classic() +
  labs(y = "Subtype %", 
       x = "Treatment Status")+
  rremove("legend")+
  theme(axis.text = element_text(size=18,angle = 0, hjust = 1),
        axis.title = element_text(size=20),
        axis.text.x = element_text(angle = 0, hjust = .5),
        strip.text = element_text(size=20))

p3

tiff("figures/subtype_treatment_status_percent_results.tiff", width=300,height=100, units = "mm", res=600)
print(p3)
dev.off()

################################################################################
# Tarla Status
################################################################################
data_df <- sce@colData %>% 
  as.data.frame() %>% 
  select(subtype,tarla,sample_id) %>% 
  filter(!is.na(tarla))

wide_data_df <- data_df %>% 
  count(subtype,sample_id) %>% 
  group_by(sample_id) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n/total)*100) %>% 
  select(subtype,sample_id,freq) %>% 
  pivot_wider(names_from = "subtype",values_from = "freq",values_fill = 0) %>% 
  merge(.,data_df,by='sample_id') %>% 
  select(tarla,A,N,P,I,sample_id) %>% 
  distinct()


wide_data_df$tarla <- factor(wide_data_df$tarla, levels=c("pre","post"))

curr_subtype <- "N"

results_list <- list()
pvals <- c()
for(curr_subtype in c("A","N","P","I")){
  
  formula_str <- glue("{curr_subtype} ~ {colnames(wide_data_df)[1]}")
  
  model <- lm(formula_str, data = wide_data_df)
  
  model_summary <- summary(model)
  
  curr_pval <- model_summary$coefficients[2,4]
  
  pvals <- append(pvals,curr_pval)
  
  results_list <- append(results_list, list(model))
  
}


names(pvals) <- c("A","N","P","I")

pval_df <- pvals %>% 
  enframe() %>% 
  rename("subtype"="name","pval"="value")

plot_df <- wide_data_df %>% 
  pivot_longer(!c(sample_id,tarla), names_to = "subtype",values_to = "freq") %>% 
  merge(.,pval_df,by="subtype")

plot_df$subtype <- factor(plot_df$subtype, levels = c("A","N","P","I"))
plot_df$tarla <- ifelse(plot_df$tarla == "pre","Pre-Tarlatamab","Post-Tarlatamab")

p4 <- ggplot(plot_df, aes(x = tarla, y = freq)) +
  facet_wrap(~subtype,nrow=1)+
  geom_boxplot(aes(fill=tarla),outlier.shape = NA) +
  geom_text(aes(label=glue("p-value = {sprintf('%.2f',pval)}"),x=1.5,y=103,size=5))+
  geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
  scale_fill_manual(values=c("#A8DADC", "#dd4b33"))+
  theme_classic() +
  labs(y = "Subtype %", 
       x = "Treatment Status")+
  rremove("legend")+
  theme(axis.text = element_text(size=18,angle = 0, hjust = 1),
        axis.title = element_text(size=20),
        axis.text.x = element_text(angle = 0, hjust = .5,size=12),
        strip.text = element_text(size=20))

p4

tiff("figures/subtype_tarla_percent_results.tiff", width=350,height=100, units = "mm", res=600)
print(p4)
dev.off()
