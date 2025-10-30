library(glue)
library(reshape2)
library(readxl)
library(HDCytoData)
library(CATALYST)
library(ComplexHeatmap)
library(ggpubr)
library(diffcyt)
library(lme4)
library(tidytext)
library(RColorBrewer)
library(circlize)
library(limma)
library(factoextra)
library(iMUBAC)
library(cyCombine)
library(egg)
library(gridExtra)
library(rstatix) 
library(ggalluvial)
library(broom)
library(broom.mixed)
library(purrr)

library(tidyverse)

scale_values <- function(x){(x-min(x))/(max(x)-min(x))}

cytof_de <- function(sce, method = "wilcox", metric = "median", ident = "new_clusters"){
  
  state_markers <- as.data.frame(rowData(sce)) %>% 
    dplyr::filter(marker_class == "state") %>% 
    pull(marker_name)
  
  # counts <- sce@assays@data$counts
  counts <- sce@assays@data$exprs
  
  counts <- counts[state_markers,]
  
  idents <- as.character(unique(colData(sce)[[ident]]))
  
  protein <- c()
  ident_list <- c()
  logfc <- c()
  p_val <- c()
  
  if(method == "wilcox"){
    
    for(i in 1:length(state_markers)){
      
      for(curr_ident in idents){
        protein <- append(protein, state_markers[i])
        ident_list <- append(ident_list, curr_ident)
        
        a <- counts[i,colData(sce)[[ident]] == curr_ident]
        b <- counts[i,colData(sce)[[ident]] != curr_ident]
        
        wilcox_res <- wilcox.test(a,b)
        
        p_val <- append(p_val,wilcox_res$p.value)
        
        pseudocount <- 1e-6 
        
        #calculate log2FC
        if(metric == "mean"){
          logfc <- append(logfc, log2((mean(a)+pseudocount)/(mean(b)+pseudocount)))
        } else if(metric == "median"){
          logfc <- append(logfc, log2((median(a)+pseudocount)/(median(b)+pseudocount)))
        }
        
      }
    }
  }
  
  if(method == "glm"){
    
    for(i in 1:length(state_markers)){
      
      # Set up data for fitting model
      data_i <- cbind(counts[i,],colData(sce))
      colnames(data_i)[1] <- "protein"
      
      # Fit glm
      fit <- glm(glue("protein ~ {ident} + experiment_id + patient_id + condition"), data = data_i)
      
      # Extract p values for idents
      cluster_pvals <- summary(fit)$coefficients[,4][1:length(idents)]
      p_val <- append(p_val,cluster_pvals)
      
      for(curr_ident in idents){
        protein <- append(protein, state_markers[i])
        ident_list <- append(ident_list, curr_ident)
        
        #calculate log2FC
        a <- counts[i,colData(sce)[[ident]] == curr_ident]
        b <- counts[i,colData(sce)[[ident]] != curr_ident]
        
        if(metric == "mean"){
          logfc <- append(logfc, log2(mean(a)/mean(b)))
        } else if(metric == "median"){
          logfc <- append(logfc, log2(median(a)/median(b)))
        }
        
      }
      
    }
  }
  
  df <- data.frame(cbind(protein,ident_list,logfc,p_val))
  
  df$p_adj <- p.adjust(df$p_val, method = "BH")
  
  df$logfc <- as.numeric(df$logfc)
  
  #remove NaNs and Inf
  df <- df %>%
    dplyr::filter(!is.infinite(logfc) & !is.nan(logfc))
  
  df <- df %>%
    dplyr::filter(abs(logfc) > .05)
  
  return(df)
}

create_expression_heatmap <- function(sce, group, markers, curr_title="", scale=T){
  
  all_groups <- sort(unique(sce[[group]]))
  
  if(group == "subtype"){
    all_groups <- c("A","N","P","Mes")
  }
  
  
  heatmap <- matrix(NA, ncol=length(markers),nrow=length(all_groups))
  for(i in 1:length(all_groups)){
    
    curr_group <- all_groups[i]
    # cat(curr_group,"\n")
    
    for(j in 1:length(markers)){
      
      curr_marker <- markers[j]
      # cat(curr_marker,"\n")
      
      heatmap[i,j] <- mean(sce@assays@data$exprs[curr_marker,sce[[group]] == curr_group])
    }
  }
  
  colnames(heatmap) <- markers
  rownames(heatmap) <- all_groups
  
  
  if(scale){
    # scaled_heatmap <- t(scale(t(heatmap)))
    scaled_heatmap <- scale(heatmap)
  } else{
    scaled_heatmap <- heatmap
  }
  
  
  return(scaled_heatmap)
 

  
  # return(t(scaled_heatmap))
}

create_marker_boxplots <- function(sce, markers_to_use, group, fill = NULL, alpha = NULL, nrow_to_use=3){
  
  y <- assay(sce, "exprs")
  
  df <- data.frame(t(y), colData(sce), check.names = FALSE)
  
  value <- ifelse("exprs" == "exprs", "expression", "exprs")
  
  gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
                id.vars = names(colData(sce)))
  
  ################################################################################
  
  if(group == "subtype"){
    gg_df$subtype <- factor(gg_df$subtype, levels=c("A","N","P",'Mes'))
  }
  
  plot_df <- gg_df %>%
    dplyr::filter(antigen %in% markers_to_use)
  
  plot_df$antigen <- factor(plot_df$antigen, levels = markers_to_use)
  
  if(is.null(fill)){
    p <- ggboxplot(plot_df, x=group ,y="expression", fill=group, outlier.size = .1)
    # p <- ggviolin(plot_df, x=group ,y="expression", fill=group, outlier.size = .1)
  } else{
    p <- ggboxplot(plot_df, x=group ,y="expression", fill=fill, outlier.size = .1)
    # p <- ggviolin(plot_df, x=group ,y="expression", fill=group, outlier.size = .1)
  }
  
  if(group == "cancer_enriched"){
    p <- p+theme(axis.text.x = element_text(size=10),
                 legend.position = "top")
  } 
  
  
  p1 <- p + facet_wrap(~antigen,nrow=nrow_to_use)+
    theme(axis.title = element_text(size=16),
      strip.text = element_text(face = "bold", size=14), 
           strip.background = element_blank())
  
  return(p1)
  
}

# diff_abundance_barpl <- ot <- function(sce, group1, group2){
#   
#   all_group1 <- unique(sce[[group1]])
#   all_group2 <- unique(sce[[group2]])
#   
#   pvals <- c()
#   ORs <- c()
#   
#   for(curr_group1 in all_group1){
#     
#     a <- sum(colData(treatment_status_sce)$new_clusters == curr_group1 & colData(treatment_status_sce)$treatment_status == "treated")
#     b <- sum(colData(treatment_status_sce)$new_clusters == curr_group1 & colData(treatment_status_sce)$treatment_status != "treated")
#     c <- sum(colData(treatment_status_sce)$new_clusters != curr_group1 & colData(treatment_status_sce)$treatment_status == "treated")
#     d <- sum(colData(treatment_status_sce)$new_clusters != curr_group1 & colData(treatment_status_sce)$treatment_status != "treated")
#     
#     contin_table <- matrix(c(a+.5,c+.5,b+.5,d+.5),ncol=2)
#     
#     fisher_res <- fisher.test(contin_table)
#     
#     ORs <- append(ORs,fisher_res$estimate)
#     pvals <- append(pvals,fisher_res$p.value)
#   }
#   
#   # Select significant clusters
#   signif_clusters <- which(p.adjust(pvals) < 0.05)
#   
#   cluster_prop_df <- as.data.frame(colData(treatment_status_sce)) %>% 
#     dplyr::count(new_clusters,treatment_status) %>% 
#     group_by(treatment_status) %>% 
#     mutate(total = sum(n)) %>% 
#     mutate(freq = (n / total)*100)
#   
#   # Add star for significance
#   cluster_prop_df <- cluster_prop_df %>% 
#     mutate(significant = ifelse(new_clusters %in% signif_clusters, "*","")) %>% 
#     group_by(new_clusters) %>% 
#     mutate(height = max(freq))
#   
#   
#   cluster_prop_df$tarla <- ifelse(cluster_prop_df$treatment_status == "treated", "Treated", "Naive")
#   
#   condition_colors <- c("Treated" = "#FFB74D","Naive"="#64B5F6")
#   
#   p1 <- ggplot(cluster_prop_df,aes(x=new_clusters,y=freq,fill=tarla))+
#     geom_col(position = "dodge")+
#     geom_text(aes(y = height+.01,label=significant),size=6)+
#     xlab("Cluster")+
#     ylab("Percentage")+
#     labs(fill="Condition")+
#     scale_fill_manual(values=condition_colors)+
#     theme_classic()+
#     theme(panel.grid.minor = element_blank(), 
#           strip.text = element_text(face = "bold", size=8), 
#           axis.text = element_text(color = "black", size=8),
#           axis.title = element_text(size=8),
#           legend.text = element_text(size=6),
#           legend.title = element_text(size=8))
#   
#   
#   return(p1)
# }

# Write sessionInfo to text file
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")