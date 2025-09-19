source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################
ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

################################################################################
# Treatment status and subtype association
################################################################################

contin_table <- table(ctcs$sex, ctcs$subtype) %>% 
  as.data.frame.matrix()

# colnames(contin_table) <- c("ASCL1","Inflamed","NeuroD1","POU2F3")


# contin_table <- contin_table[c(1,3),]
# contin_table <- contin_table[c(2,4),]

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

res_df$subtype <- factor(res_df$subtype, levels=c("I","P","N",'A'))

status_plot <- ggplot(res_df, aes(x = status, y = subtype)) +
  geom_point(aes(size = abs(residual), fill = residual), shape = 21, color = "black", stroke = 0.6) +
  geom_text(aes(label = stars), vjust = -1.5, size = 4) +  # Stars above dots
  scale_fill_gradient2(low = "#457B9D", mid = "white", high = "#dd4b33", midpoint = 0) +
  scale_size(range = c(2, 8)) +
  theme_minimal()+
  theme(axis.text = element_text(size=12,angle = 0, hjust = 1),
        axis.title = element_text(size=14),
        axis.text.x = element_text(angle = 0, hjust = .5))+
  labs(title = "",
       size = "|Residual|", fill = "Residual",
       y = "Subtype",
       x = "")

status_plot
