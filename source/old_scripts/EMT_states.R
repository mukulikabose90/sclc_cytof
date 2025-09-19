source("source/cytof_de_function.R")

script_seed <- 42

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

ctcs <- ctcs[rowData(ctcs)$marker_class == "state",]

EMT_proteins <- c("E-Cad","SLUG","Vimentin","Twist","EpCAM","MUC-1","CD24",'CD44')

rownames(ctcs)

ctcs <- ctcs[EMT_proteins,]

ctcs <- CATALYST::cluster(ctcs, features = "state",
                          xdim = 10, ydim = 10, maxK = 20, seed = script_seed)


ctcs <- runDR(ctcs, "UMAP", cells = 5e3, features = "state")


ctcs@metadata$delta_area

plotDR(ctcs, "UMAP", color_by = "meta9", facet_by = "patient_id",scale = T)+
  geom_point(size=1)



plotDR(ctcs, "UMAP", color_by = "treatment_status",scale = T)+
  geom_point(size=1)

colData(ctcs)$new_clusters <- cluster_ids(ctcs, "meta9")

# Plot UMAP manually
xy <- reducedDim(ctcs, "UMAP")[, 1:2]
colnames(xy) <- c("x", "y")
df <- data.frame(colData(ctcs), xy, check.names = FALSE)

display.brewer.pal(5,"Spectral")

# Generate and save cluster colors
cluster_colors <- colorRampPalette(brewer.pal(9,"Spectral"))(10)

# Plot UMAP
p1 <- ggplot(df)+
  geom_point(aes(x=x, y=y, color=new_clusters), size= 0.5)+
  xlab("UMAP 1")+
  ylab("UMAP 2")+
  scale_color_manual(name = "Clusters",values=cluster_colors)+
  theme_classic() +
  guides(color = guide_legend(override.aes = list(size=5)))+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=20), 
        axis.text = element_text(color = "black", size=15),
        axis.title = element_text(size=20),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18))

p1

png("figures/ctcs_emt_states_umap.png", width = 8, height = 8, units = "in", res = 1200)
print(p1)
dev.off()

ggplot(df)+
  geom_point(aes(x=x, y=y, color=subtype))+
  xlab("UMAP 1")+
  ylab("UMAP 2")+
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size=5)))+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=20), 
        axis.text = element_text(color = "black", size=15),
        axis.title = element_text(size=20),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18))

################################################################################
# Cluster expression boxplots
y <- assay(ctcs, "exprs")

df <- data.frame(t(y), colData(ctcs), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")
gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(ctcs)))

p2 <- ggplot(gg_df)+
  geom_boxplot(aes(x=new_clusters, y=expression,fill=new_clusters))+
  facet_wrap(~antigen, scales="free_y",nrow=2)+
  scale_fill_manual(values=cluster_colors)+
  guides(fill="none")+
  xlab("Clusters")+
  ylab("Expression")+
  theme_classic()+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=20), 
        axis.text = element_text(color = "black", size=15),
        axis.title = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=20))

png("figures/ctcs_emt_states_boxplots.png", width = 12, height = 6, units = "in", res = 1200)
print(p2)
dev.off()


ggplot(gg_df)+
  geom_boxplot(aes(x=antigen, y=expression,fill=antigen))+
  facet_wrap(~new_clusters, scales="free_y",nrow=3)+
  scale_fill_manual(values=cluster_colors)+
  guides(fill="none")+
  xlab("Clusters")+
  ylab("Expression")+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=20), 
        axis.text = element_text(color = "black", size=15),
        axis.text.x = element_text(color = "black", size=15, angle=45,hjust=1,vjust=1),
        axis.title = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=20))


################################################################################
# Subtype proportions in each cluster
cluster_prop_df <- as.data.frame(colData(ctcs)) %>%
  group_by(new_clusters,subtype) %>%
  summarise(n = n()) %>%
  group_by(new_clusters) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n / total*100))

cluster_prop_df <- cluster_prop_df %>% 
  dplyr::filter(total > 5)

ggplot(cluster_prop_df,aes(x=new_clusters, y=freq, fill=subtype))+
  geom_col()+
  geom_text(aes(y = 103,label=total))+
  xlab("Cluster")+
  ylab("Percentage")+
  labs(fill="Subtype")+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=20), 
        axis.text = element_text(color = "black", size=15),
        axis.title = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=20))

################################################################################
# cluster proportions in each subtype
cluster_prop_df <- as.data.frame(colData(ctcs)) %>%
  group_by(new_clusters,subtype) %>%
  summarise(n = n()) %>%
  group_by(subtype) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n / total)*100)

cluster_prop_df <- cluster_prop_df %>% 
  dplyr::filter(total > 5)

ggplot(cluster_prop_df,aes(x=subtype, y=freq, fill=new_clusters))+
  geom_col(position = "dodge")+
  geom_text(aes(y = 103,label=total))+
  xlab("Subtype")+
  ylab("Percentage")+
  labs(fill="Cluster")+
  scale_fill_manual(values=cluster_colors)+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=20), 
        axis.text = element_text(color = "black", size=15),
        axis.title = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=20))

################################################################################
# Subtype proportions in each cluster
cluster_prop_df <- as.data.frame(colData(ctcs)) %>%
  group_by(new_clusters,treatment_status) %>%
  summarise(n = n()) %>%
  group_by(new_clusters) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n / total*100))

cluster_prop_df <- cluster_prop_df %>% 
  dplyr::filter(total > 5)

# Capitalize status
cluster_prop_df$treatment_status <- ifelse(cluster_prop_df$treatment_status == "naive", "Naive","Treated")


ggplot(cluster_prop_df,aes(x=new_clusters, y=freq, fill=treatment_status))+
  geom_col()+
  geom_text(aes(y = 103,label=total))+
  xlab("Cluster")+
  ylab("Percentage")+
  scale_fill_manual(values=c("deepskyblue3", "darkorange"))+
  labs(fill="Treatment Status")+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=20), 
        axis.text = element_text(color = "black", size=15),
        axis.title = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=20))
#############################################
# Cluster DE

df <- cytof_de(ctcs, method = "wilcox", metric = "median", ident = "new_clusters")

plot_df <- df %>% 
  mutate(significance = factor(ifelse(p_adj < .05, "signif", "ns"), levels=c("signif","ns"))) %>% 
  mutate(protein = reorder_within(protein, logfc, ident_list))

plot_df <- df %>% 
  mutate(significance = factor(ifelse(p_adj < .05, "*", ""), levels=c("*",""))) %>% 
  mutate(protein = reorder_within(protein, logfc, ident_list)) %>% 
  mutate(star_x = ifelse(logfc > 0, logfc+.25,logfc-.5)) %>% 
  mutate(up_down= ifelse(logfc>0,"up","down"))


# Plot DE results of significant clusters
ggplot(plot_df,aes(x=as.numeric(logfc), y=protein, fill=up_down))+
  geom_col(color="darkgray",size=.001)+
  geom_text(aes(x=star_x, label=significance), size = 5)+
  facet_wrap(~ident_list, scales = "free_y")+
  scale_y_reordered()+
  labs(fill = "")+
  guides(fill="none")+
  scale_fill_manual(values=c("royalblue", "firebrick1"))+
  ylab("Protein")+
  xlab("Median log(FC)")+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=20), 
        axis.text = element_text(color = "black", size=15),
        axis.title = element_text(size=20),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18))


##########################################################################################

clusters <- levels(colData(ctcs)$new_clusters)

pvals <- c()
ORs <- c()

for(curr_cluster in clusters){
  
  a <- sum(colData(ctcs)$new_clusters == curr_cluster & colData(ctcs)$treatment_status == "naive")
  b <- sum(colData(ctcs)$new_clusters == curr_cluster & colData(ctcs)$treatment_status != "naive")
  c <- sum(colData(ctcs)$new_clusters != curr_cluster & colData(ctcs)$treatment_status == "naive")
  d <- sum(colData(ctcs)$new_clusters != curr_cluster & colData(ctcs)$treatment_status != "naive")
  
  contin_table <- matrix(c(a+.5,c+.5,b+.5,d+.5),ncol=2)
  
  fisher_res <- fisher.test(contin_table)
  
  ORs <- append(ORs,fisher_res$estimate)
  pvals <- append(pvals,fisher_res$p.value)
}

# Select significant clusters
signif_clusters <- which(p.adjust(pvals) < 0.05)

cluster_prop_df <- as.data.frame(colData(ctcs)) %>% 
  dplyr::filter(treatment_status %in% c("naive","treated")) %>% 
  dplyr::count(new_clusters,treatment_status) %>% 
  group_by(treatment_status) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n / total)*100)

# Add star for significance
cluster_prop_df <- cluster_prop_df %>% 
  mutate(significant = ifelse(new_clusters %in% signif_clusters, "*","")) %>% 
  group_by(new_clusters) %>% 
  mutate(height = max(freq))


# Capitalize status
cluster_prop_df$treatment_status <- ifelse(cluster_prop_df$treatment_status == "naive", "Naive","Treated")


p3 <- ggplot(cluster_prop_df,aes(x=new_clusters,y=freq,fill=treatment_status))+
  geom_col(position = "dodge")+
  geom_text(aes(y = height+.01,label=significant),size=10)+
  xlab("Cluster")+
  ylab("Percentage")+
  scale_fill_manual(values=c("deepskyblue3", "darkorange"))+
  labs(fill="Treatment Status")+
  theme_classic()+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=20), 
        axis.text = element_text(color = "black", size=15),
        axis.title = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=20))

png("figures/ctcs_emt_states_treatment_stage_barplots.png", width = 12, height = 6, units = "in", res = 1200)
print(p3)
dev.off()

##########################################################################################





