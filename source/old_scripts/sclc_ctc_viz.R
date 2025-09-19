# This script is the original script I created to test visualization methods using
# CyTOF data.

library(tidyverse)
library(glue)
library(reshape2)
library(readxl)
library(HDCytoData)
library(CATALYST)
library(ComplexHeatmap)
library(ggpubr)

##################################################################################

sce <- readRDS("data/cytof_objects/sclc_cytof_sce_object.rds")

sce <- sce[,colData(sce)$condition != "cell_line"]
# PLOT EXPRESSION DENSITY PLOTS

ggplot(gg_df, fill = NULL, aes(x = value, y = "after_stat(ndensity)", col = patient_id, group = "sample_id")) + 
  facet_wrap(~antigen, scales = "free_x") +
  geom_density() + 
  ylab("normalized density") + 
  theme_classic() +
  theme(panel.grid = element_blank(), 
                          strip.background = element_blank(), strip.text = element_text(face = "bold"), 
                          axis.text = element_text(color = "black"), axis.title = element_text(color = "black"))


Heatmap(y)

p <- plotExprs(sce, color_by = "patient_id")
p$facet$params$ncol <- 6
p

# PLOT MEDIAN EXPRESSION HEATMAP

plotExprHeatmap(sce, scale = "last",by='patient_id')

y <- assay(sce, "exprs")
df <- data.frame(t(y), colData(sce), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")
gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(sce)))

median_heatmap <- gg_df %>% 
  select(patient_id,antigen,expression) %>% 
  group_by(patient_id,antigen) %>% 
  summarize(median_value = median(expression)) %>% 
  pivot_wider(names_from = antigen,values_from = median_value) %>% 
  column_to_rownames("patient_id")

Heatmap(median_heatmap)


# PLOT CELL COUNTS BAR PLOT
plotCounts(sce, group_by="patient_id")+
  xlab("Patient")+
  ylab("Cell Count")


# PLOT PCA

pbMDS(sce, color_by = "patient_id", label_by = "patient_id")


          

# plotNRS(sce, features = "none", color_by = "condition")


y <- assay(sce, "exprs")

# Scale data. Scale protein expression between 0-1 across all cells
scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
y <- t(apply(y, MARGIN = 1, FUN = function(X) scale_values(X)))

df <- data.frame(t(y), colData(sce), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")
gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(sce)))

ggplot(gg_df)+
  geom_boxplot(aes(x=condition, y=expression))+
  facet_wrap(~antigen, scales="free_y")+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))


ggplot(gg_df)+
  geom_boxplot(aes(x=patient_id, y=expression))+
  facet_wrap(~antigen, scales="free_y")+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))




unique(gg_df$patient_id)

p <- ggboxplot(gg_df %>% 
                 dplyr::filter(antigen == "POU2F3"), x = "condition", y = "expression",
               color = "condition", palette = "jco",
               add = "jitter")


p + stat_compare_means(method = "wilcox.test")

unique(sce$experiment_id)

my_comparisons <- list( c("513549", "514521"), c("513549", "515600"), c("515600", "514521") )
ggboxplot(gg_df %>% 
            dplyr::filter(patient_id == "SC414"), x = "experiment_id", y = "expression",
          color = "condition", palette = "jco",
          add = "jitter") + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons)


corr_df <- gg_df %>% 
  dplyr::filter(antigen %in% c("DLL3", "ASCL1") & condition == "cancer") %>%
  pivot_wider(names_from = antigen,
              values_from = expression) %>% 
  na.omit()


str(corr_df)
sum(is.na(corr_df))

ggplot(corr_df,aes(y=ASCL1,x=DLL3))+
  geom_point() + 
  facet_wrap(~patient_id)+
  geom_smooth(method = "lm", se = T)

cor.test(corr_df$DLL3,corr_df$ASCL1, method = "pearson")





