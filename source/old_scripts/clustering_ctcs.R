source("source/cytof_de_function.R")

script_seed <- 42

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

ctcs <- ctcs[rowData(ctcs)$marker_class == "state",]

ctcs <- CATALYST::cluster(ctcs, features = "state",
                          xdim = 10, ydim = 10, maxK = 20, seed = script_seed)

plotDR(ctcs, "UMAP", color_by = "meta4",scale = T)+
  geom_point(size=1)

plotDR(ctcs, "UMAP", color_by = "subtype",scale = T)+
  geom_point(size=1)

plotDR(ctcs, "UMAP", color_by = c("NeuroD1","ASCL1","POU2F3"),scale = T)+
  geom_point(size=1)


###########################################

EMT_proteins <- c("NeuroD1","ASCL1","POU2F3")

rownames(ctcs)

ctcs <- ctcs[EMT_proteins,]

ctcs <- CATALYST::cluster(ctcs, features = "state",
                          xdim = 10, ydim = 10, maxK = 20, seed = script_seed)


ctcs <- runDR(ctcs, "UMAP", cells = 5e3, features = "state")


ctcs@metadata$delta_area

plotDR(ctcs, "UMAP", color_by = "subtype",scale = T)+
  geom_point(size=1)

plotDR(ctcs, "UMAP", color_by = "patient_id",scale = T)+
  geom_point(size=1)



all_patients <- as.character(unique(ctcs$patient_id))

# Generate and save patient colors
patient_colors <- colorRampPalette(brewer.pal(length(all_patients),"Spectral"))(length(all_patients))
names(patient_colors) <- all_patients

saveRDS(patient_colors, "data/patient_colors.rds")

colData(ctcs)
xy <- reducedDim(ctcs, "UMAP")[, 1:2]
colnames(xy) <- c("x", "y")
df <- data.frame(colData(ctcs), xy, check.names = FALSE)

# Plot UMAP
ggplot(df)+
  geom_point(aes(x=x, y=y, color=treatment_status))+
  xlab("UMAP 1")+
  ylab("UMAP 2")+
  scale_color_manual(name = "Patient ID",values=patient_colors)+
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size=5)))+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=20), 
        axis.text = element_text(color = "black", size=15),
        axis.title = element_text(size=20),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18))

curr_patient <- all_patients[1]
for(curr_patient in all_patients){
  
  curr_patient_colors <- ifelse(names(patient_colors) == curr_patient, patient_colors,"gray")
  
  curr_patient_alpha <- ifelse(names(patient_colors) == curr_patient, 1,.05)
  
  xy <- reducedDim(ctcs, "UMAP")[, 1:2]
  colnames(xy) <- c("x", "y")
  df <- data.frame(colData(ctcs), xy, check.names = FALSE)
  
  ggplot(df)+
    geom_point(aes(x=x, y=y, color=patient_id, alpha = patient_id))+
    xlab("UMAP 1")+
    ylab("UMAP 2")+
    scale_alpha_manual(name = "Patient ID",values=curr_patient_alpha)+
    scale_color_manual(name = "Patient ID",values=curr_patient_colors)+
    theme_minimal() +
    guides(color = guide_legend(override.aes = list(size=5)))+
    theme(panel.grid.minor = element_blank(), 
          strip.text = element_text(face = "bold", size=20), 
          axis.text = element_text(color = "black", size=15),
          axis.title = element_text(size=20),
          legend.text = element_text(size=15),
          legend.title = element_text(size=18))
}


ctcs$tre

plotDR(ctcs, "UMAP", color_by = c("NeuroD1","ASCL1","POU2F3"),scale = T)+
  geom_point(size=1)


