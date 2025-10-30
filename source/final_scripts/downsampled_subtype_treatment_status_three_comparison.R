source("source/sclc_cytof_functions.R")

set.seed(42)
################################################################################
# Read in data
################################################################################
ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

cluster_colors <- c("#dd4b33", "#D1DACF", "#A8DADC", "#457B9D")

################################################################################
# Sample n cells from each patient
################################################################################
n_cells <- 30

sampled_data <- ctcs@colData %>% 
  as.data.frame() %>%
  group_by(collection_id) %>%
  filter(n() >= n_cells) %>%        # keep only patients with ≥ n_cells cells
  slice_sample(n = n_cells) %>%
  ungroup()

sampled_data$treatment_status <- ifelse(sampled_data$treatment_status == "naive","Naive","CTX ± ICI")

sampled_data$treatment_status <- ifelse(is.na(sampled_data$tarla), sampled_data$treatment_status,
                                        ifelse(sampled_data$tarla == "pre", sampled_data$treatment_status, "Tarla"))

################################################################################
treatment_table <- sampled_data %>% 
  count(subtype,treatment_status) %>% 
  arrange(subtype) %>% 
  pivot_wider(names_from = treatment_status,values_from = n) %>% 
  column_to_rownames("subtype") 

subtypes <- rownames(treatment_table)

results_list <- list()
i <- 1
for(i in 1:4){
  
  contin_table <- rbind(treatment_table[i,],colSums(treatment_table[-i,]))
  
  # Chi-squared test: for significance
  chi_result <- chisq.test(contin_table, correct = F)
  chi_result
  
  ################################################################################
  # Chi-squared test: Naive vs CTX ± ICI
  tbl <- contin_table[,c(2,1)]
  chi_result <- chisq.test(tbl, correct = F)
  
  a <- tbl[1,1]
  b <- tbl[1,2]
  c <- tbl[2,1]
  d <- tbl[2,2]
  
  OR <- (a*d)/(b*c)
  
  # Standard error of log(OR)
  SE_logOR <- sqrt(1/a + 1/b + 1/c + 1/d)
  
  # 95% CI
  CI_lower <- exp(log(OR) - 1.96*SE_logOR)
  CI_upper <- exp(log(OR) + 1.96*SE_logOR)
  
  results_list <- append(results_list, list(c(subtypes[i],"Naive_CTX ± ICI",OR,CI_upper,CI_lower,chi_result$p.value)))
  
  ################################################################################
  # Chi-squared test: Naive vs Tarla
  tbl <- contin_table[,c(3,1)]
  chi_result <- chisq.test(tbl, correct = F)
  
  a <- tbl[1,1]
  b <- tbl[1,2]
  c <- tbl[2,1]
  d <- tbl[2,2]
  
  OR <- (a*d)/(b*c)
  
  # Standard error of log(OR)
  SE_logOR <- sqrt(1/a + 1/b + 1/c + 1/d)
  
  # 95% CI
  CI_lower <- exp(log(OR) - 1.96*SE_logOR)
  CI_upper <- exp(log(OR) + 1.96*SE_logOR)

  results_list <- append(results_list, list(c(subtypes[i],"Naive_Tarla",OR,CI_upper,CI_lower,chi_result$p.value)))
  
  ################################################################################
  # Chi-squared test: CTX ± ICI vs Tarla
  tbl <- contin_table[,c(3,2)]
  chi_result <- chisq.test(tbl, correct = F)
  
  a <- tbl[1,1]
  b <- tbl[1,2]
  c <- tbl[2,1]
  d <- tbl[2,2]
  
  OR <- (a*d)/(b*c)
  
  # Standard error of log(OR)
  SE_logOR <- sqrt(1/a + 1/b + 1/c + 1/d)
  
  # 95% CI
  CI_lower <- exp(log(OR) - 1.96*SE_logOR)
  CI_upper <- exp(log(OR) + 1.96*SE_logOR)
  
  ################################################################################
  results_list <- append(results_list, list(c(subtypes[i],"CTX ± ICI_Tarla",OR,CI_upper,CI_lower,chi_result$p.value)))

}

plot_df <- as.data.frame(do.call(rbind,results_list))
colnames(plot_df) <- c("subtype","comparison","or","upper_or", "lower_or", "pval")
plot_df$or <- as.numeric(plot_df$or)
plot_df$upper_or <- as.numeric(plot_df$upper_or)
plot_df$lower_or <- as.numeric(plot_df$lower_or)
plot_df$pval <- as.numeric(plot_df$pval)

################################################################################

plot_df <- plot_df %>% 
  mutate(padj = p.adjust(pval), method = "BH",) %>% 
  mutate(signif = ifelse(padj < 0.05, "s","ns")) %>% 
  mutate(log_or = log(or)) %>% 
  mutate(log_upper_or = log(upper_or)) %>% 
  mutate(log_lower_or = log(lower_or))

plot_df$left_label <- as.character(sapply(plot_df$comparison, function(x) strsplit(x,"_")[[1]][1]))
plot_df$right_label <- as.character(sapply(plot_df$comparison, function(x) strsplit(x,"_")[[1]][2]))

plot_df$left_label[plot_df$left_label=="CTX ± ICI"] <- "CTX ± ICI"
plot_df$right_label[plot_df$right_label=="CTX ± ICI"] <- "CTX ± ICI"


plot_df$subtype <- factor(plot_df$subtype, levels=c("A","N","P","Mes"))
plot_df$comparison <- factor(plot_df$comparison, levels = c("Naive_CTX ± ICI","Naive_Tarla","CTX ± ICI_Tarla"))
################################################################################

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



tiff(glue("figures/downsampled_subtype_treatment_status_chi_squared_results_all_comp_{n_cells}.tiff"), width=360,height=200, units = "mm", res=600)
print(p1)
dev.off()



