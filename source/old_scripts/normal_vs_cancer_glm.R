


model_df <- data.frame(t(y),colData(sce)) %>% 
  select(-c(cell_id,sample_id,experiment_id,patient_id))

model_df$condition <- factor(model_df$condition, levels=c("normal","cancer"))


expr_model <- glm("condition ~ .", data = model_df, family = binomial(link = logit))

temp <- as.data.frame(summary(expr_model)$coefficients)

colnames(temp) <- c("estimate", "stderr","z","pval")

temp <- temp %>% 
  rownames_to_column("antigen") %>% 
  dplyr::filter(antigen != "(Intercept)") %>% 
  dplyr::filter(pval < 0.05) %>% 
  arrange(desc(estimate))

temp$antigen <- factor(temp$antigen, levels=temp$antigen)
ggplot(temp)+
  geom_col(aes(x=antigen,y=estimate, fill=estimate))+
  scale_fill_gradient2(low="blue", mid= "white",high="red")


