source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################

files_to_read <- list.files("data/run_stats/")

all_data <- list()

for(i in 1:length(files_to_read)){
  data <- read.csv(glue("data/run_stats/{files_to_read[i]}"))
  
  
  data <- data[,c(which(grepl("FCS",colnames(data))),which(grepl("pop",colnames(data))),which(grepl("experiment_id",colnames(data))))]
  colnames(data) <- c("filename","event_count","experiment_id")
  
  all_data <- append(all_data,list(data))
  
}

all_data <- do.call(rbind,all_data)

# Remove cell lines
all_data <- all_data[!(grepl("NJH",all_data$filename)|
             grepl("H1105",all_data$filename)|
             grepl("H526",all_data$filename)|
             grepl("MRC5",all_data$filename)|
               grepl("H196",all_data$filename)|
               grepl("Unassign",all_data$filename)|
               grepl("PE",all_data$filename)),]

# Subset to only cancer samples
all_data <- all_data[grepl("SC",all_data$filename),]

all_data$filename <- sapply(all_data$filename, FUN = function(x) strsplit(x, "_")[[1]][1])

all_data$filename[which(all_data$filename == "SC293")] <- "SC293-4"

all_data <- all_data %>% 
  mutate(sample_id = paste0(filename,"_",experiment_id))


# Read in ctc data
ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

ctcs_count <- ctcs %>% 
  colData() %>% 
  as.data.frame() %>% 
  count(sample_id)


all_data <- merge(all_data,ctcs_count,by="sample_id")

all_data <- all_data %>% 
  mutate(percent_ctc = sprintf("%.3f",(n/event_count)*100))


all_data <- merge(all_data,as.data.frame(colData(ctcs)),by="sample_id")


  

summary(all_data$pct_ctc)

ctc_table <- all_data %>% 
  select(sample_id,percent_ctc) %>% 
  distinct() 


write.csv(ctc_table,file = "data/ctc_percent_table.csv",row.names = F)





plot_df <- all_data %>% 
  select(pct_ctc,treatment_status) %>% 
  distinct()




vec1 <- plot_df %>% 
  filter(treatment_status == "naive") %>% 
  pull(pct_ctc)
  

vec2 <- plot_df %>% 
  filter(treatment_status == "treated") %>% 
  pull(pct_ctc)


boxplot(vec1,vec2)
