################################################################################
# Read in data
################################################################################

sce <- readRDS("data/cytof_objects/sclc_all_samples_object_no_qc.rds")

################################################################################
# Get protein markers
################################################################################
marker_info <- read.csv("data/cytof_panel_info.csv")
marker_info <- data.frame(marker_info, stringsAsFactors = FALSE)

state_markers <- marker_info %>%
  dplyr::filter(marker_class == "state") %>%
  pull(antigen)

saveRDS(state_markers, "data/state_markers.rds")

################################################################################
# Remove outlier experiment
################################################################################
sce <- sce[,sce$experiment_id != "531050"]

sce <- sce[,sce$sample_type != "pleural_effusion"]


y <- assay(sce, "exprs")

df <- data.frame(t(y), colData(sce), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")

gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(sce)))


markers_to_use <- c("p-Rb","DLL3","MUC-1","Alcam")
markers_to_use <- c("Ir191","Ir193")

plot_df <- gg_df %>%
  dplyr::filter(antigen %in% markers_to_use)


ggviolin(plot_df, x="condition" ,y="expression", fill="sample_type")+
  facet_grid(antigen~experiment_id)



all_experiments <- unique(sce$experiment_id)

all_experiments <- sce@colData %>% 
  as.data.frame() %>% 
  select(experiment_id,sample_type, condition) %>% 
  distinct() %>% 
  count(experiment_id) %>% 
  filter(n > 2) %>% 
  pull(experiment_id) %>% 
  as.character()



curr_experiment <- all_experiments[1]

plot_df <- gg_df %>% 
  filter(experiment_id == curr_experiment & antigen %in% markers_to_use[1])


plot_df <- gg_df %>% 
  filter(experiment_id %in% all_experiments & antigen %in% markers_to_use&sample_type!="cell_line") %>% 
  mutate(group = paste0(sample_type,"_",condition))


ggplot(plot_df)+
  geom_density(aes(x=expression, fill=group), alpha=.3)+
  facet_grid(antigen~experiment_id)



curr_experiment <- all_experiments[1]

plot_df <- gg_df %>% 
  filter(experiment_id == curr_experiment & antigen %in% markers_to_use[1]) %>% 
  mutate(group = paste0(sample_type,"_",condition))



cancer_expr <- plot_df %>% 
  filter(group == "blood_cancer") %>% 
  pull(expression)

normal_expr <- plot_df %>% 
  filter(group == "blood_normal") %>% 
  pull(expression)

cell_line_expr <- plot_df %>% 
  filter(group == "cell_line_cancer") %>% 
  pull(expression)



get_cv(cell_line_expr)
get_cv(cancer_expr)
get_cv(normal_expr)

max_val <- mean(normal_expr)+(sd(normal_expr)*2)
min_val <- mean(normal_expr)-(sd(normal_expr)*2)

sum(cancer_expr > max_val | cancer_expr < min_val)


get_cv <- function(data){
  data_mean <- mean(data)
  data_sd <- sd(data)
  
  cv <- (data_sd / data_mean) * 100
  
  return(cv)
}




