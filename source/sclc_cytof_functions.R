library(tidyverse)
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
        
        #calculate log2FC
        if(metric == "mean"){
          logfc <- append(logfc, log2(mean(a)/mean(b)))
        } else if(metric == "median"){
          logfc <- append(logfc, log2(median(a)/median(b)))
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
    dplyr::filter(abs(logfc) > .25)
  
  return(df)
}
