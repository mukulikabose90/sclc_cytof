source("source/cytof_de_function.R")

script_seed <- 42

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

ctcs <- ctcs[rowData(ctcs)$marker_class == "state",]

EMT_proteins <- c("NeuroD1","ASCL1","POU2F3")

ctcs <- ctcs[EMT_proteins,]

ctcs <- CATALYST::cluster(ctcs, features = "state",
                          xdim = 10, ydim = 10, maxK = 20, seed = script_seed)


ctcs <- runDR(ctcs, "UMAP", cells = 5e3, features = "state")


all_patients <- as.character(unique(ctcs$patient_id))


as.data.frame(colData(ctcs)) %>% 
  dplyr::filter(treatment_status == "") %>% 
  pull(patient_id) %>% 
  unique()

# Generate and save patient colors
patient_colors <- colorRampPalette(brewer.pal(length(all_patients),"Spectral"))(length(all_patients))
names(patient_colors) <- all_patients


colData(ctcs)
xy <- reducedDim(ctcs, "UMAP")[, 1:2]
colnames(xy) <- c("x", "y")
df <- data.frame(colData(ctcs), xy, check.names = FALSE)




# Plot UMAP
p1 <- ggplot(df)+
  geom_point(aes(x=x, y=y, color=treatment_status))+
  xlab("UMAP 1")+
  ylab("UMAP 2")+
  scale_color_manual(name = "Treatment Status",values=c("gray","royalblue","orangered"),labels = c("Unknown","Naive","Treated"))+
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size=5)))+
  theme(panel.grid.minor = element_blank(), 
        strip.text = element_text(face = "bold", size=20), 
        axis.text = element_text(color = "black", size=15),
        axis.title = element_text(size=20),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18))

dir.create("figures/patient_ctc_umaps/", showWarnings = FALSE)
png("figures/patient_ctc_umaps/treatment_stage_umap.png", width = 8000,height = 6000, res = 1000)
p1
dev.off()
