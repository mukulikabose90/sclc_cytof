source("source/cytof_de_function.R")

sce <- readRDS("data/cytof_objects/sclc_cytof_sce_object_with_clusters.rds")

sce <- sce[,colData(sce)$condition == "cancer"]

sce <- sce[rowData(sce)$marker_class == "state",]

sce <- sce[,colData(sce)$new_clusters %in% c(1,2)]


y <- assay(sce, "exprs")

# Scale data. Scale protein expression between 0-1 across all cells
# scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
# y <- t(apply(y, MARGIN = 1, FUN = function(X) scale_values(X)))

df <- data.frame(t(y), colData(sce), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")
gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(sce)))

ggplot(gg_df)+
  geom_boxplot(aes(x=treatment_status, y=expression))+
  facet_wrap(~antigen, scales="free_y")+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))


p <- ggboxplot(gg_df, x = "treatment_status", y = "expression",
               color = "treatment_status", palette = "jco", facet.by = "antigen")


p + stat_compare_means(method = "wilcox.test",label.y = 7.5)

################################


df <- cytof_de(sce, method = "wilcox", metric = "median", ident = "treatment_status")


plot_df <- df %>% 
  mutate(significance = factor(ifelse(p_adj < .05, "signif", "ns"), levels=c("signif","ns"))) %>% 
  mutate(protein = reorder_within(protein, logfc, ident_list))


# Plot DE results of significant clusters
ggplot(plot_df) +
  geom_col(aes(x=as.numeric(logfc), y=protein, fill=significance))+
  facet_wrap(~ident_list, scales = "free_y")+
  scale_y_reordered()+
  ylab("Protein")+
  xlab("Mean log2FC")




state_markers <- as.data.frame(rowData(sce)) %>% 
  dplyr::filter(marker_class == "state") %>% 
  pull(marker_name)

counts <- sce@assays@data$exprs

counts <- counts[state_markers,]

clusters <- c("naive","treated")
protein <- c()
cluster <- c()
logfc <- c()
p_val <- c()
metric <- "median"
for(i in 1:length(state_markers)){
  
  for(curr_cluster in clusters){
    protein <- append(protein, state_markers[i])
    cluster <- append(cluster, curr_cluster)
    
    a <- counts[i,colData(sce)$treatment_status == curr_cluster]
    b <- counts[i,colData(sce)$treatment_status != curr_cluster]
    
    wilcox_res <- wilcox.test(a,b)
    
    p_val <- append(p_val,wilcox_res$p.value)
    
    #calculate log2FC
    if(metric == "mean"){
      logfc <- append(logfc, log2(mean(a)/mean(b)))
    } else if(metric == "median"){
      logfc <- append(logfc, log2(median(a)/median(b)))
    }
    
  }
}

df <- data.frame(cbind(protein,cluster,logfc,p_val))

df$p_adj <- p.adjust(df$p_val, method = "BH")

df$logfc <- as.numeric(df$logfc)

#remove NaNs and Inf
df <- df %>% 
  dplyr::filter(!is.infinite(logfc) & !is.nan(logfc))

df <- df %>%
  dplyr::filter(abs(logfc) > .25)


plot_df <- df %>% 
  mutate(significance = factor(ifelse(p_adj < .05, "signif", "ns"), levels=c("signif","ns"))) %>% 
  mutate(protein = reorder_within(protein, logfc, cluster))

ggplot(plot_df)+
  geom_col(aes(x=as.numeric(logfc), y=protein, fill=significance))+
  facet_wrap(~cluster, scales = "free_y")+
  scale_y_reordered()+
  ylab("Protein")+
  xlab("Median log2FC")


