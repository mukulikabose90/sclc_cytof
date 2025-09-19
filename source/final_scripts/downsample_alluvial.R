################################################################################
# This script downsamples cells from patients to n_cells and then plots the 
# proportions of subtypes and treatment status as an alluvial plot.
################################################################################
source("source/sclc_cytof_functions.R")

set.seed(42)
################################################################################
# Read in data
################################################################################
ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

cluster_colors <- c("#dd4b33", "#D1DACF", "#A8DADC", "#457B9D")

################################################################################
# # Sample n cells from each patient
################################################################################
n_cells <- 20
for(n_cells in c(20,50)){
  sampled_data <- ctcs@colData %>% 
    as.data.frame() %>%
    group_by(patient_id) %>%
    filter(n() >= n_cells) %>%        # keep only patients with ≥ n_cells cells
    slice_sample(n = n_cells) %>%
    ungroup()
  
  cat(length(unique(sampled_data$patient_id)), " patients\n")
  
  ################################################################################
  # Create plot dataframe
  ################################################################################
  plot_df <- sampled_data %>% 
    select(treatment_status,subtype,tarla) %>% 
    cbind(1) %>% 
    rename("n" = "1")
  
  plot_df$treatment_status <- ifelse(plot_df$treatment_status == "naive","Naive","SOC")
  plot_df$tarla <- ifelse(plot_df$tarla == "pre","Pre-Tarla","Tarla")
  
  plot_df$treatment_status <- ifelse(plot_df$treatment_status == "Naive", plot_df$treatment_status, ifelse(plot_df$tarla == "Pre-Tarla" | is.na(plot_df$tarla),plot_df$treatment_status,"Tarla"))
  
  plot_df$treatment_status <- factor(plot_df$treatment_status, levels=c("Naive","SOC","Pre-Tarla","Tarla"))
  
  plot_df$subtype <- factor(plot_df$subtype, levels=c("A","N","P",'I'))
  
  plot_df_long <- to_lodes_form(data.frame(plot_df),
                                key = "category", value = "group", id = "cohort",
                                axes = 1:2) %>% 
    add_count(group) %>% 
    group_by(category) %>% 
    mutate(total=sum(n)) %>% 
    mutate(freq=sprintf("%.1f",(nn/total)*100)) 
  
  
  plot_df_long_left <- subset(plot_df_long, group %in% c("Naive","SOC","Tarla"))
  plot_df_long_right <- subset(plot_df_long, !group %in% c("Naive","SOC","Tarla"))
  
  ################################################################################
  # Plot alluvial 
  ################################################################################
  p1 <- ggplot(data = plot_df_long,
               aes(x = category, stratum = group, alluvium = cohort, y = total)) +
    coord_flip() +
    scale_y_reverse() +
    geom_flow(aes(fill=group),width=.3, aes.flow = "backward") +
    geom_stratum(aes(fill=group),width=.3) +
    geom_text(stat = "stratum", aes(label = glue("{group}")),size=12) +
    ggrepel::geom_text_repel(data=plot_df_long_left,stat = "stratum", aes(label = glue("n={nn}\n({freq}%)")),nudge_x = -.3,
                             size=10,segment.color = NA) +
    ggrepel::geom_text_repel(data=plot_df_long_right,stat = "stratum", aes(label = glue("n={nn}\n({freq}%)")),nudge_x = .3,
                             size=10,segment.color = NA) +
    scale_fill_manual(name = "group",values=c("gray90","gray8 0","gray70",cluster_colors))+
    theme_void() +
    rremove("legend")
  
  ################################################################################
  # Treatment Status Association
  ################################################################################
  data_df <- sampled_data
  
  data_df$treatment_status <- factor(data_df$treatment_status, levels = c("naive","treated"))
  
  results_list <- list()
  for(curr_subtype in c("A","N","P","I")){
    data_df$curr_subtype <- as.factor(as.integer(data_df$subtype == curr_subtype))
    
    formula_str <- glue("curr_subtype ~ treatment_status + (1 | patient_id)")
    
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
    mutate(padj = p.adjust(pval)) %>% 
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
  
  ################################################################################
  # Tarla Status Association
  ################################################################################
  data_df <- sampled_data
  
  data_df$tarla <- factor(data_df$tarla, levels = c("pre","post"))
  
  results_list <- list()
  for(curr_subtype in c("A","N","P","I")){
    data_df$curr_subtype <- as.factor(as.integer(data_df$subtype == curr_subtype))
    
    formula_str <- glue("curr_subtype ~ tarla + (1 | patient_id)")
    
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
    mutate(padj = p.adjust(pval)) %>% 
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
    xlim(-2.5,2.5)+
    labs(y="Subtype",
         x="log(OR)")+
    theme_classic()+
    annotate("text", x=-.75, y=4.5, label = "Pre-Tarlatamab", angle=0,size=4) +
    annotate("text", x=.75, y=4.5, label = "Post-Tarlatamab", angle=0,size=4) +
    theme(axis.text = element_text(size=22,angle = 0, hjust = 1),
          axis.title = element_text(size=24),
          axis.text.x = element_text(angle = 0, hjust = .5)) 
  
  ################################################################################
  # Save figures
  ################################################################################
  
  tiff(glue("figures/cell_level_alluvial_plot_downsampled_{n_cells}.tiff"), width=500,height=300, units = "mm", res=600)
  print(p1)
  dev.off()
  
  tiff(glue("figures/downsampled_subtype_treatment_status_or_results_{n_cells}.tiff"), width=200,height=200, units = "mm", res=600)
  print(p3)
  dev.off()
  
  tiff(glue("figures/downsampled_subtype_tarla_or_results_{n_cells}.tiff"), width=200, height=200, units = "mm", res=600)
  print(p4)
  dev.off()
}







