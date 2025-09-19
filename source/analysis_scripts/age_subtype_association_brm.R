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

# age_df$age[age_df$age < 40] <- 1
# age_df$age[age_df$age > 39 & age_df$age < 50] <- 2
# age_df$age[age_df$age > 49 & age_df$age < 60] <- 3
# age_df$age[age_df$age > 59 & age_df$age < 70] <- 4 
# age_df$age[age_df$age > 69 & age_df$age < 80] <- 5
# age_df$age[age_df$age > 79] <- 6
# 
# age_df$age <- factor(age_df$age, levels=c(1,2,3,4,5,6))

age_df$subtype <- factor(age_df$subtype, levels=c("A","N","P",'I'))
age_df$sex <- factor(age_df$sex)


################################################################################
# Test for significance
model <- brm(
  formula = subtype ~ age + (1 | sample_id),
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


################################################################################
# Visualize

ce <- conditional_effects(model, effects = "age", categorical = TRUE)
t <- plot(ce)

df <- ce$age

head(df,10)

ggplot(df, aes(x = age, y = estimate__, color = effect2__, fill = effect2__)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, color = NA) +
  scale_color_manual(values = c("#dd4b33", "#F1FAEE", "#A8DADC", "#457B9D"))+
  scale_fill_manual(values = c("#dd4b33", "#F1FAEE", "#A8DADC", "#457B9D"))+
  labs(
    x = "Age",
    y = "Predicted Probability",
    color = "Subtype",
    fill = "Subtype",
    title = "Predicted Probability of Subtype by Age") +
    theme_classic()
