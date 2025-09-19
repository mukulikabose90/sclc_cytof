

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

normals <- readRDS("data/cytof_objects/normals_with_subtype.rds")

final_df <- list()
for(i in c("ASCL1","NeuroD1","POU2F3")){
  temp <- ctcs[,ctcs$subtype == substring(i, 1,1)]
  temp2 <- normals[,normals$subtype == substring(i, 1,1)]
  
  
  df <- as.data.frame(rbind(cbind(temp@assays@data$counts[i,],"ctc"),
                            cbind(temp2@assays@data$counts[i,],"normal")))
  
  df <- cbind(df,i)
  
  colnames(df) <- c("value","class","protein")
  df$value <- as.numeric(df$value)
  
  
  final_df <- append(final_df, list(df))
  
}


final_df <- do.call(rbind,final_df)

ggplot(df)+
  geom_density(aes(x=value,fill=class))+
  facet_wrap(~protein)

ggplot(final_df)+
  geom_boxplot(aes(y=value,fill=class))+
  facet_wrap(~protein, scales="free_y")

