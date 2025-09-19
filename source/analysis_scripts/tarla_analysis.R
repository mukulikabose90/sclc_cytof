source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################

sce <- readRDS("data/cytof_objects/all_samples_ctcs_with_subtype.rds")

colData(sce)$condition <- factor(colData(sce)$condition, levels=c("normal", "cancer"))
sce@metadata$experiment_info$condition <- factor(sce@metadata$experiment_info$condition, levels=c("normal", "cancer"))

sce <- sce[,colData(sce)$condition == "cancer"]

sce$tarla <- as.character(sce$tarla)

sce <- sce[,(!is.na(sce$tarla)) & sce$tarla != "NA "]

tarla_sce <- sce



dim(tarla_sce)
table(tarla_sce$tarla)

tarla_sce <- CATALYST::cluster(tarla_sce, features = "state",
                         xdim = 10, ydim = 10, maxK = 20, seed = script_seed)


tarla_sce <- runDR(tarla_sce, "UMAP", cells = 5e3, features = "state")

CATALYST::plotDR(tarla_sce, color_by = "meta4")

tarla_sce@metadata$delta_area

colData(tarla_sce)$new_clusters <- cluster_ids(tarla_sce, "meta4")

saveRDS(tarla_sce, "data/cytof_objects/tarla_sce.rds")

# Plot UMAP manually
xy <- reducedDim(tarla_sce, "UMAP")[, 1:2]
colnames(xy) <- c("x", "y")
df <- data.frame(colData(tarla_sce), xy, check.names = FALSE)


# cluster_colors <- c(
#   "#E57373",  # muted red
#   "#64B5F6",  # muted blue
#   "#81C784",  # muted green
#   "#FFB74D"   # muted orange
# )

cluster_colors <- c(
  "#E57373",  # muted red
  "#FFB74D",  # muted orange
  "#BA68C8",  # muted purple
  "#64B5F6",  # muted blue
  "#81C784"  # muted green
)


# Plot UMAP
p1 <- ggplot(df)+
  geom_point(aes(x=x, y=y, color=new_clusters),size=3)+
  xlab("UMAP 1")+
  ylab("UMAP 2")+
  labs(color = "Clusters")+
  # scale_color_manual(name = "Clusters",values=cluster_colors)+
  theme_classic() +
  guides(color = guide_legend(override.aes = list(size=5)))+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=8), 
        axis.text = element_text(color = "black", size=8),
        axis.title = element_text(size=8),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8))



p1

df$tarla <- factor(df$tarla, levels=c("pre","post"))

facet_names <- c('pre'="Pre-Tarlatamab",'post'="Post-Tarlatamab")
p2 <- ggplot(df)+
  geom_point(aes(x=x, y=y, color=new_clusters),size=3)+
  facet_wrap(~tarla,labeller=as_labeller(facet_names))+
  xlab("UMAP 1")+
  ylab("UMAP 2")+
  labs(color = "Clusters")+
  # scale_color_manual(name = "Clusters",values=cluster_colors)+
  guides(color = guide_legend(override.aes = list(size=5)))+
  scale_alpha_manual(values = c("ctc" = 1, "non-ctc" = 0.05))+
  theme_classic() +
  guides(alpha = "none")+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=8), 
        axis.text = element_text(color = "black", size=8),
        axis.title = element_text(size=8),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8))  

p2

jpeg("figures/tarla_cluster.jpg", width=120,height=100, units = "mm", res=1000)
print(p1)
dev.off()


jpeg("figures/tarla_pre_vs_post_cluster.jpg", width=180,height=100, units = "mm", res=1000)
print(p2)
dev.off()


################################################################################

clusters <- levels(tarla_sce$new_clusters)

pvals <- c()
ORs <- c()

for(curr_cluster in clusters){
  
  a <- sum(colData(tarla_sce)$new_clusters == curr_cluster & colData(tarla_sce)$tarla == "T")
  b <- sum(colData(tarla_sce)$new_clusters == curr_cluster & colData(tarla_sce)$tarla != "T")
  c <- sum(colData(tarla_sce)$new_clusters != curr_cluster & colData(tarla_sce)$tarla == "T")
  d <- sum(colData(tarla_sce)$new_clusters != curr_cluster & colData(tarla_sce)$tarla != "T")
  
  contin_table <- matrix(c(a+.5,c+.5,b+.5,d+.5),ncol=2)
  
  fisher_res <- fisher.test(contin_table)
  
  ORs <- append(ORs,fisher_res$estimate)
  pvals <- append(pvals,fisher_res$p.value)
}

# Select significant clusters
signif_clusters <- which(p.adjust(pvals) < 0.05)

cluster_prop_df <- as.data.frame(colData(tarla_sce)) %>% 
  dplyr::count(new_clusters,tarla) %>% 
  group_by(tarla) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n / total)*100)

# Add star for significance
cluster_prop_df <- cluster_prop_df %>% 
  mutate(significant = ifelse(new_clusters %in% signif_clusters, "*","")) %>% 
  group_by(new_clusters) %>% 
  mutate(height = max(freq))


cluster_prop_df$tarla <- ifelse(cluster_prop_df$tarla == "post", "Post-Tarlatamab", "Pre-Tarlatamab")

cluster_prop_df$tarla <- factor(cluster_prop_df$tarla, levels=c("Pre-Tarlatamab","Post-Tarlatamab"))

condition_colors <- c("Pre-Tarlatamab" = "#F48FB1","Post-Tarlatamab"="#BA68C8")

p1 <- ggplot(cluster_prop_df,aes(x=new_clusters,y=freq,fill=tarla))+
  geom_col(position = "dodge")+
  geom_text(aes(y = height+.01,label=significant),size=6)+
  xlab("Cluster")+
  ylab("Percentage")+
  labs(fill="Condition")+
  scale_fill_manual(values=condition_colors)+
  theme_classic()+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=8), 
        axis.text = element_text(color = "black", size=8),
        axis.title = element_text(size=8),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8))


p1

jpeg("figures/tarla_cluster_diff_abundance_barplots.jpg", width=120,height=100, units = "mm", res=1000)
print(p1)
dev.off()
################################################################################
subtypes <- unique(colData(tarla_sce)$subtype)

pvals <- c()
ORs <- c()

for(curr_subtype in subtypes){
  
  a <- sum(colData(tarla_sce)$subtype == curr_subtype & colData(tarla_sce)$tarla == "T")
  b <- sum(colData(tarla_sce)$subtype == curr_subtype & colData(tarla_sce)$tarla != "T")
  c <- sum(colData(tarla_sce)$subtype != curr_subtype & colData(tarla_sce)$tarla == "T")
  d <- sum(colData(tarla_sce)$subtype != curr_subtype & colData(tarla_sce)$tarla != "T")
  
  contin_table <- matrix(c(a+.5,c+.5,b+.5,d+.5),ncol=2)
  
  fisher_res <- fisher.test(contin_table)
  
  ORs <- append(ORs,fisher_res$estimate)
  pvals <- append(pvals,fisher_res$p.value)
}

# Select significant subtypes
signif_subtypes <- subtypes[which(p.adjust(pvals) < 0.05)]

subtype_prop_df <- as.data.frame(colData(tarla_sce)) %>% 
  dplyr::count(subtype,tarla) %>% 
  group_by(tarla) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n / total)*100)

# Add star for significance
subtype_prop_df <- subtype_prop_df %>% 
  mutate(significant = ifelse(subtype %in% signif_subtypes, "*","")) %>% 
  group_by(subtype) %>% 
  mutate(height = max(freq))


subtype_prop_df$tarla <- ifelse(subtype_prop_df$tarla == "post", "Post-Tarlatamab", "Pre-Tarlatamab")

subtype_prop_df$tarla <- factor(subtype_prop_df$tarla, levels=c("Pre-Tarlatamab","Post-Tarlatamab"))

condition_colors <- c("Pre-Tarlatamab" = "#F48FB1","Post-Tarlatamab"="#BA68C8")

p1 <- ggplot(subtype_prop_df,aes(x=subtype,y=freq,fill=tarla))+
  geom_col(position = "dodge")+
  geom_text(aes(y = height+.01,label=significant),size=6)+
  xlab("Subtype")+
  ylab("Percentage")+
  labs(fill="Condition")+
  scale_fill_manual(values=condition_colors)+
  theme_classic()+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=8), 
        axis.text = element_text(color = "black", size=8),
        axis.title = element_text(size=8),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8))


p1

jpeg("figures/tarla_subtype_diff_abundance_barplots.jpg", width=120,height=100, units = "mm", res=1000)
print(p1)
dev.off()

################################################################################
signif_clusters <- c(1,2)

metric_to_use <- "mean"
df <- cytof_de(tarla_sce, method = "wilcox", metric = metric_to_use, ident = "new_clusters")


plot_df <- df %>% 
  mutate(significance = factor(ifelse(p_adj < .05, "*", ""), levels=c("*",""))) %>% 
  mutate(protein = reorder_within(protein, logfc, ident_list)) %>% 
  mutate(star_x = ifelse(logfc > 0, logfc+.25,logfc-.5)) %>% 
  mutate(up_down= ifelse(logfc>0,"up","down"))


plot_df$ident_list <- paste0("Cluster ",plot_df$ident_list)
plot_df$ident_list <- factor(plot_df$ident_list, levels=paste0("Cluster ", c(1:10)))


signif_df <- plot_df %>% 
  dplyr::filter(ident_list %in% paste0("Cluster ", signif_clusters))


x_axis_label <- gsub("m","M",metric_to_use)

p2 <- ggplot(plot_df,aes(x=as.numeric(logfc), y=protein, fill=logfc))+
  geom_col(color="darkgray",size=.001)+
  geom_text(aes(x=star_x, label=significance), size = 3)+
  facet_wrap(~ident_list, scales = "free", nrow=2)+
  scale_y_reordered()+
  labs(fill = "")+
  guides(fill="none")+
  # scale_fill_manual(values=c("royalblue", "firebrick1"))+
  scale_fill_gradient2(low = "royalblue4", mid = "white", high = "firebrick", midpoint = 0)+
  ylab("Protein")+
  xlab(glue("{x_axis_label} log(FC)"))+
  theme_classic()+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=8), 
        axis.text = element_text(color = "black", size=6),
        axis.title = element_text(size=8),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8))

p2

jpeg(glue("figures/tarla_cluster_{metric_to_use}_de_all.jpg"), width=140,height=100, units = "mm", res=1000)
print(p2)
dev.off()





