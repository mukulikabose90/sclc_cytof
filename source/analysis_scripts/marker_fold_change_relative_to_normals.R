


sce <- readRDS("data/cytof_objects/sclc_all_samples_object.rds")

y <- assay(sce, "exprs")

df <- data.frame(t(y), colData(sce), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")

gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(sce)))




all_samples <- as.character(unique(gg_df$sample_id))


all_samples <- all_samples[!grepl("NORMAL",all_samples)]

markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam", "E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist", "SLUG", "PD-L1", "p-YAP", "CD44", "CD24")

all_markers <- markers_to_use

sample_median_fc_values <- list()

curr_sample <- all_samples[45]
for(curr_sample in all_samples){
  curr_experiment <- strsplit(curr_sample, "_")[[1]][2]
  
  curr_conditions <- gg_df %>% 
    filter(experiment_id == curr_experiment) %>% 
    count(condition) %>% 
    pull(condition)
  
  # continue if this sample was run with normals
  if("normal" %in% curr_conditions){
    
    for(curr_marker in all_markers){
      
      normal_expr_vector <- gg_df %>% 
        filter(antigen == curr_marker & condition == "normal" & experiment_id == curr_experiment) %>% 
        pull(expression)
      
      normal_median_expr <- median(normal_expr_vector)
      
      expr_vector <- gg_df %>% 
        filter(antigen == curr_marker & sample_id == curr_sample & experiment_id == curr_experiment) %>% 
        pull(expression)
      
      median_expr <- median(expr_vector)
      
    
      curr_fc_val <- median_expr+.001/normal_median_expr+.001
      
      
      sample_median_fc_values <- append(sample_median_fc_values,list(c(curr_sample,curr_marker,curr_fc_val)))
      
      
    }
  }
}

final_df <- as.data.frame(do.call(rbind,sample_median_fc_values))
colnames(final_df) <- c("sample_id","antigen","fc")
final_df$fc <- log(as.numeric(final_df$fc))



ggviolin(final_df, x="antigen",y="fc",draw_quantiles = .5)+
  geom_hline(yintercept=0)+
  labs(y="log(FC)")





