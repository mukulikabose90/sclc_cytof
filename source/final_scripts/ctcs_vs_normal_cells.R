source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################
library(ggbeeswarm)
ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")


sce <- readRDS("data/cytof_objects/sclc_all_samples_object.rds")


sce$ctc <- ifelse(sce$cell_id %in% ctcs$cell_id, "ctc","normal")


y <- assay(sce, "exprs")

df <- data.frame(t(y), colData(sce), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")

gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
               id.vars = names(colData(sce)))


markers_to_use <- c("p-Rb","DLL3","MUC-1","Alcam")
markers_to_use <- c("Ir191","Ir193")

plot_df <- gg_df %>%
  dplyr::filter(antigen %in% markers_to_use)


p <- ggviolin(plot_df, x="ctc" ,y="expression", fill="ctc")

ggboxplot(plot_df, x="ctc",y="expression",fill="ctc",add = "jitter", add.params = list(size = 1, alpha = 0.025))+
  facet_wrap(~antigen)+
  stat_compare_means(label = "p.signif", label.y = c(7,11,7,9), tip.length = 0, size=8, comparisons = list(c("normal","ctc")))+
  labs(y=glue("Expression"),
       x="")+
  ylim(c(0,12))+
  rremove("legend")+
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=20))


p1 <- p + facet_wrap(~antigen, scales = "free_y")+
  stat_compare_means()

p1




rows_to_keep <- intersect(rownames(sce),rownames(ctcs))
sce <- sce[rows_to_keep,]


ctcs$cell_id
sce$cell_id
sce$sample_id



ctcs <- ctcs[,ctcs$treatment_status == "naive"]
normal_cells <- sce[,sce$condition == "normal"]


curr_marker <- "Alcam"



plot_df <- as.data.frame(rbind(cbind(ctcs@assays@data$exprs[curr_marker,],"CTCs"),cbind(normal_cells@assays@data$exprs[curr_marker,],"Normal Cells")))


colnames(plot_df) <- c("value","class")

plot_df$value <- as.numeric(plot_df$value)

ggviolin(plot_df, x="class",y="value",fill="class",add = "none", add.params = list(size = 1, alpha = 0.3))+
  stat_compare_means(label = "p.signif", label.y = 6, tip.length = 0, size=8, comparisons = list(c("CTCs","Normal Cells")))+
  labs(y=glue("{curr_marker} Expression"),
       x="")+
  ylim(c(0,6.5))+
  rremove("legend")+
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=20))




  geom_beeswarm(binaxis = "y",
              binwidth = .05,
              stackdir = "center",
              dotsize = 0.01)




  
  
  y <- assay(ctcs, "exprs")
  
  df <- data.frame(t(y), colData(ctcs), check.names = FALSE)
  
  value <- ifelse("exprs" == "exprs", "expression", "exprs")
  
  gg_df2 <- melt(df, value.name = "expression", variable.name = "antigen", 
                 id.vars = names(colData(ctcs)))
  
  dim(gg_df1)
  dim(gg_df2)
  
  
  
  ################################################################################
  
  plot_df <- gg_df %>%
    dplyr::filter(antigen %in% markers_to_use)
  
  
  if(is.null(fill)){
    p <- ggboxplot(plot_df, x=group ,y="expression", fill=group)
  } else{
    p <- ggboxplot(plot_df, x=group ,y="expression", fill=fill)
  }
  
  
  p1 <- p + facet_wrap(~antigen, scales = "free_y")+
    stat_compare_means()
  


