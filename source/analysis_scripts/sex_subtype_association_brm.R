source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)

cluster_colors <- c("#dd4b33", "#F1FAEE", "#A8DADC", "#457B9D")
################################################################################
ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

age_df <- ctcs@colData %>% 
  as.data.frame() %>% 
  filter(!is.na(age) & !is.na(sex)) %>% 
  select(subtype,age,sex,sample_id)

age_df$age <- as.numeric(as.character(age_df$age))

age_df$subtype <- factor(age_df$subtype, levels=c("A","N","P",'I'))
age_df$sex <- factor(age_df$sex)

plot_df <- age_df %>% 
  filter(sex %in% c("Male","Female")) %>% 
  dplyr::count(sex,subtype) %>% 
  group_by(sex) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n/total)*100) 

plot_df$subtype <- factor(plot_df$subtype, levels=c("A","N","P",'I'))
plot_df$sex <- factor(plot_df$sex, levels=c("Male","Female"))

ggplot(plot_df)+
  geom_col(aes(x=sex,y=freq,fill=subtype))+
  scale_fill_manual(name = "Subtype",values=cluster_colors)+
  geom_text(aes(label=total,x=sex), y=105,size = 3)+
  ylim(0,105)+
  scale_fill_manual(name = "Subtype",values=cluster_colors)+
  labs(x="",
       y="Percentage",
       fill="Subtype")+
  theme_classic()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=10, color="black"),
        strip.text = element_text(face = "bold", size=12), 
        strip.background = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12))

################################################################################
# Test for significance
model <- brm(
  formula = subtype ~ sex + (1 | sample_id),
  family = categorical(),
  data = age_df,
  iter = 100,
  cores = 12,
  chains = 4)

fit_null <- brm(
  formula = subtype ~ 1 + (1 | sample_id), 
  family = categorical(), 
  data = age_df,
  iter = 100,
  cores = 12,
  chains = 4)

loo_full <- loo(model)
loo_null <- loo(fit_null)

loo_compare(loo_full, loo_null)
