


ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

plot_df <- ctcs@colData %>% 
  as.data.frame() 

# plot_df$sample_id <- sapply(as.character(plot_df$sample_id), function(x) strsplit(x, "_")[[1]][1])

sample_level_subtype <- plot_df %>% 
  select(treatment_status, subtype, collection_id) %>% 
  group_by(collection_id) %>% 
  count(subtype) %>% 
  filter(n == max(n)) %>% 
  select(collection_id,subtype) 

tied_samples <- sample_level_subtype %>% 
  count(collection_id) %>% 
  filter(n > 1)

df <- plot_df %>% 
  select(collection_id,treatment_status) %>% 
  merge(.,sample_level_subtype, by="collection_id") %>% 
  distinct() %>% 
  select(collection_id,treatment_status,subtype)


plot_df <- as.data.frame(table(df)) %>% 
  filter(treatment_status != "normal") %>% 
  group_by(treatment_status) %>% 
  mutate(total = sum(Freq)) %>% 
  mutate(freq = (Freq/total)*100)

plot_df$subtype <- factor(plot_df$subtype, levels = c("A","N","P","I"))

ggplot(plot_df)+
  geom_col(aes(x=treatment_status,y=freq,fill=subtype))+
  geom_text(aes(x=treatment_status,y=101,label=total),size=3)+
  scale_fill_manual(values = cluster_colors)+
  theme_classic()



contin_table <- as.data.frame.matrix(table(df))
contin_table <- contin_table[c(1,3),]

fisher.test(contin_table)


test_res <- chisq.test(contin_table)

residuals <- test_res$stdres

# Calculate two-sided p-values from standard normal distribution
pvals <- 2 * (1 - pnorm(abs(residuals)))

# Adjust using Benjamini-Hochberg (FDR)
pvals_adj <- matrix(p.adjust(as.vector(pvals), method = "BH"),
                    nrow = nrow(pvals),
                    dimnames = dimnames(pvals))

# Convert to long format
res_df <- as.data.frame(as.table(residuals)) %>%
  mutate(p_value = as.vector(pvals_adj),
         stars = case_when(
           p_value < 0.001 ~ "***",
           p_value < 0.01  ~ "**",
           p_value < 0.05  ~ "*",
           TRUE            ~ ""))

colnames(res_df) <- c("status","subtype","residual","adj_pval","stars")
res_df$status <- ifelse(res_df$status == "post", "Post-Tarlatamab","Pre-Tarlatamab")

# res_df$subtype <- factor(res_df$subtype, levels=c("A","N","P",'I'))
res_df$subtype <- factor(res_df$subtype, levels=c("I","P","N",'A'))

res_df$status <- factor(res_df$status, levels=c("Pre-Tarlatamab","Post-Tarlatamab"))

tarla_plot <- ggplot(res_df, aes(x = status, y = subtype)) +
  geom_point(aes(size = abs(residual), fill = residual), shape = 21, color = "black", stroke = 0.6) +
  geom_text(aes(label = stars), vjust = -1.5, size = 4) +  # Stars above dots
  scale_fill_gradient2(low = "#457B9D", mid = "white", high = "#dd4b33", midpoint = 0) +
  scale_size(range = c(2, 8)) +
  theme_minimal() +
  theme(axis.text = element_text(size=12,angle = 0, hjust = 1),
        axis.title = element_text(size=14),
        axis.text.x = element_text(size=10,angle = 0, hjust = .5))+
  labs(title = "",
       size = "|Residual|", fill = "Residual",
       y = "Subtype",
       x = "")

tarla_plot

