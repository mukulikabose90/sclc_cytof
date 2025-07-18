

source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################
ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")



data_df <- ctcs@colData %>% 
  as.data.frame() %>% 
  select(subtype,age)


data_df$age <- as.numeric(as.character(data_df$age))

ggplot(data_df)+
  geom_violin(aes(x=subtype,y=age,fill=subtype))+
  geom_boxplot(aes(x=subtype,y=age,fill=subtype),width=0.1)


anova(data_df)


anova_result <- aov(age ~ subtype, data = data_df)
summary(anova_result)


ggviolin(data_df,x="subtype", y= "age",fill="subtype",draw_quantiles = 0.5)+
  stat_compare_means()
