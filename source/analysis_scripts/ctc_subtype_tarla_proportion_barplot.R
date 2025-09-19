source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################
cluster_colors <- c("#dd4b33", "#F1FAEE", "#A8DADC", "#457B9D")

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")


plot_df <- ctcs@colData %>% 
  as.data.frame() 

plot_df$sample_id <- sapply(as.character(plot_df$sample_id), function(x) strsplit(x, "_")[[1]][1])

plot_df <- plot_df %>% 
  count(sample_id,subtype,tarla) %>% 
  group_by(sample_id) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n/total)*100) %>% 
  group_by(sample_id) %>% 
  filter(!is.na(tarla))


plot_df <- plot_df %>% 
  filter(total > 10)

sample_order <- plot_df %>% 
  filter(subtype == "I") %>% 
  arrange(desc(freq)) %>% 
  pull(sample_id)

plot_df$sample_id <- factor(plot_df$sample_id, levels = sample_order)

plot_df$subtype <- factor(plot_df$subtype, levels = c("A","N","P","I"))

plot_df$tarla <- ifelse(plot_df$tarla == "pre", "Pre-Tarlatamab","Post-Tarlatamab")
plot_df$tarla <- factor(plot_df$tarla, levels = c("Pre-Tarlatamab","Post-Tarlatamab"))

ggplot(plot_df)+
  geom_col(aes(x=sample_id,y=freq,fill=subtype))+
  geom_text(aes(label=total,x=sample_id,y=101.5),size=2)+
  # facet_wrap(~tarla,scales="free")+
  facet_grid(. ~ tarla, scales = "free", space = "free") +
  scale_fill_manual(values = cluster_colors)+
  theme_classic()+
  labs(x="Sample",
       y="Percentage",
       fill="Subtype")+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))

