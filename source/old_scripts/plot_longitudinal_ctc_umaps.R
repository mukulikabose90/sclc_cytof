source("source/cytof_de_function.R")

ExpandColors <- function(colors, n, steps = 11){
  if(n <= steps){
    suppressWarnings({
      sapply(colors, function(x){colorRampPalette(c(x, "#000000"))(steps)}) %>% 
        as.data.frame() %>% 
        dplyr::filter(row_number() <= n) %>% 
        gather(key = original.color, value = expanded.color)
    })
  }else{
    warning("Select n < steps!")
  }
}

script_seed <- 42

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

ctcs <- ctcs[rowData(ctcs)$marker_class == "state",]

EMT_proteins <- c("NeuroD1","ASCL1","POU2F3")

ctcs <- ctcs[EMT_proteins,]

ctcs <- CATALYST::cluster(ctcs, features = "state",
                          xdim = 10, ydim = 10, maxK = 20, seed = script_seed)


ctcs <- runDR(ctcs, "UMAP", cells = 5e3, features = "state")

all_patients <- as.character(unique(ctcs$patient_id))

# Generate and save patient colors
patient_colors <- colorRampPalette(brewer.pal(length(all_patients),"Spectral"))(length(all_patients))
names(patient_colors) <- all_patients

# saveRDS(patient_colors, "data/patient_colors.rds")


# Get patients with longitudinal samples
long_patients <- as.data.frame(colData(ctcs)) %>% 
  dplyr::count(patient_id, collection_id) %>% 
  select(patient_id) %>% 
  table() %>% 
  as.data.frame() %>% 
  dplyr::filter(Freq > 1) %>% 
  pull(patient_id) %>% 
  as.character()

all_collections <- unique(ctcs$collection_id) %>% 
  as.character()
  

patient_collection_map <- as.data.frame(colData(ctcs)) %>% 
  select(patient_id,collection_id) %>% 
  dplyr::filter(patient_id %in% all_patients)


curr_patient <- long_patients[1]

dir.create("figures/longitudinal_ctc_umaps/", showWarnings = FALSE)
for(curr_patient in long_patients){
  
  curr_collections <- patient_collection_map %>% 
    dplyr::filter(patient_id == curr_patient) %>% 
    pull(collection_id) %>% 
    unique() %>% 
    as.character() %>% 
    sort()
  
  
  curr_colors <- ExpandColors(patient_colors[curr_patient], length(curr_collections), steps = 4)
  
  
  # show_col(curr_colors$expanded.color)
  
  curr_collection_colors <- ifelse(all_collections %in% curr_collections, rev(curr_colors$expanded.color), "gray")
  curr_collection_alpha <- ifelse(all_collections %in% curr_collections, 1, 0.05)
  
  xy <- reducedDim(ctcs, "UMAP")[, 1:2]
  colnames(xy) <- c("x", "y")
  df <- data.frame(colData(ctcs), xy, check.names = FALSE)
  
  p1 <- ggplot(df)+
    geom_point(df[df$patient_id != curr_patient,], mapping = aes(x=x, y=y), color='gray', size=0.5)+
    geom_point(df[df$patient_id == curr_patient,], mapping = aes(x=x, y=y, color=collection_id), size = 0.5)+
    xlab("UMAP 1")+
    ylab("UMAP 2")+
    ggtitle(curr_patient)+
    scale_alpha_manual(name = "Collection ID",values=curr_collection_alpha)+
    scale_color_manual(name = "Collection ID",values=(curr_colors$expanded.color))+
    theme_minimal() +
    guides(color = guide_legend(override.aes = list(size=5)))+
    theme(panel.grid.minor = element_blank(), 
          plot.title = element_text(face= "bold",size=30),
          strip.text = element_text(face = "bold", size=20), 
          axis.text = element_text(color = "black", size=15),
          axis.title = element_text(size=20),
          legend.text = element_text(size=15),
          legend.title = element_text(size=18))
  
  
  # p1
  
  
  png(paste0("figures/longitudinal_ctc_umaps/",curr_patient, "_umap.png"), width = 10000,height = 6000, res = 1000)
  print(p1)
  dev.off()
}

