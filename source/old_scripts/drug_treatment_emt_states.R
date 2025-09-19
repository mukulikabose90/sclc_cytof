source("source/cytof_de_function.R")

script_seed <- 42

ctcs <- readRDS("data/cytof_objects/drug_treatment_experiment_sce_object.rds")

ctcs <- ctcs[rowData(ctcs)$marker_class == "state",]

EMT_proteins <- c("E-Cad","SLUG","Vimentin","Twist","EpCAM","MUC-1","CD24",'CD44')

rownames(ctcs)

ctcs <- ctcs[EMT_proteins,]

ctcs <- CATALYST::cluster(ctcs, features = "state",
                          xdim = 10, ydim = 10, maxK = 20, seed = script_seed)


ctcs <- runDR(ctcs, "UMAP", cells = 5e3, features = "state")


ctcs@metadata$delta_area

plotDR(ctcs, "UMAP", color_by = "meta10",scale = T)+
  geom_point(size=1)

plotDR(ctcs, "UMAP", color_by = "meta10", facet_by="day",scale = T)+
  geom_point(size=1)


colData(ctcs)$new_clusters <- cluster_ids(ctcs, "meta10")

################################################################################
y <- assay(ctcs, "exprs")


# y <- t(scale(t(y)))
dim(y)

df <- data.frame(t(y), colData(ctcs), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")
gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(ctcs)))


gg_df$day <- factor(gg_df$day, levels=c(0,3,7,10,14))

# gg_df$new_clusters <- factor(gg_df$new_clusters, levels=c(1,2,3,4))

ggplot(gg_df)+
  geom_boxplot(aes(x=new_clusters, y=expression,fill=new_clusters))+
  facet_wrap(~antigen, scales="free_y",nrow=2)+
  guides(fill="none")+
  xlab("Clusters")+
  ylab("Expression")+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=20), 
        axis.text = element_text(color = "black", size=15),
        axis.title = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=20))



ggplot(gg_df)+
  geom_boxplot(aes(x=day, y=expression,fill=day))+
  facet_wrap(~antigen, scales="free_y",nrow=2)+
  guides(fill="none")+
  xlab("Clusters")+
  ylab("Expression")+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=20), 
        axis.text = element_text(color = "black", size=15),
        axis.title = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=20))

ggplot(gg_df)+
  geom_boxplot(aes(x=antigen, y=expression,fill=antigen))+
  facet_wrap(~new_clusters, scales="free_y",nrow=2)+
  guides(fill="none")+
  xlab("Protein")+
  ylab("Expression")+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=20), 
        axis.text = element_text(color = "black", size=15),
        axis.text.x = element_text(color = "black", size=15,angle=45,hjust=1,vjust=1),
        axis.title = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=20))

################################################################################

cluster_prop_df <- as.data.frame(colData(ctcs)) %>%
  group_by(new_clusters,day) %>%
  summarise(n = n()) %>%
  group_by(day) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n / total*100))

cluster_prop_df <- cluster_prop_df %>% 
  dplyr::filter(total > 5)

ggplot(cluster_prop_df,aes(x=new_clusters, y=freq, fill=day))+
  geom_col(position = "dodge")+
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


temp <- as.data.frame(colData(ctcs))


temp <- as.data.frame(colData(ctcs)) %>% 
  dplyr::filter(day == 0) 


temp1 <- temp[sample(rownames(temp), 248, replace = F),]

temp <- as.data.frame(colData(ctcs)) %>% 
  dplyr::filter(day == 3) 

temp2 <- temp[sample(rownames(temp), 248, replace = F),]

temp <- as.data.frame(colData(ctcs)) %>% 
  dplyr::filter(day == 7) 


temp3 <- temp[sample(rownames(temp), 248, replace = F),]

temp <- as.data.frame(colData(ctcs)) %>% 
  dplyr::filter(day == 10) 


temp4 <- temp[sample(rownames(temp), 248, replace = F),]

temp <- as.data.frame(colData(ctcs)) %>% 
  dplyr::filter(day == 14)  


temp5 <- temp[sample(rownames(temp), 248, replace = F),]


final_temp <- rbind(temp1,temp2,temp3,temp4,temp5)


rownames(temp)

cluster_prop_df <- final_temp %>%
  group_by(new_clusters,day) %>%
  summarise(n = n()) %>%
  group_by(day) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n / total*100))

cluster_prop_df <- cluster_prop_df %>% 
  dplyr::filter(total > 5)




cluster_prop_df$day <- factor(cluster_prop_df$day, levels=c(0,3,7,10,14))

ggplot(cluster_prop_df,aes(x=day, y=freq, fill=new_clusters))+
  geom_col()+
  geom_text(aes(y = 103,label=total))+
  xlab("Day")+
  ylab("Percentage")+
  labs(fill="State")+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=20), 
        axis.text = element_text(color = "black", size=15),
        axis.title = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=20))

ggplot(cluster_prop_df,aes(x=day, y=freq, color=new_clusters))+
  geom_line(aes(group=new_clusters), linewidth=1)+
  xlab("Cluster")+
  ylab("Percentage")+
  labs(fill="Subtype")+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=20), 
        axis.text = element_text(color = "black", size=15),
        axis.title = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=20))


temp <- gg_df %>% 
  group_by(day,antigen) %>% 
  summarise(med_expr = median(expression))

temp$day <- factor(temp$day, levels=c(0,3,7,10,14))

ggplot(temp)+
  geom_line(aes(x=day, y=med_expr,color=antigen, group=1))+
  facet_wrap(~antigen, scales="free_y",nrow=2)
  