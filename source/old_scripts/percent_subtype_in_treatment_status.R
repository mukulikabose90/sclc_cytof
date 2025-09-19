

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

df <- ctcs@colData %>%
  as.data.frame()

df$treatment_status <- ifelse(df$treatment_status == "naive","Naive","SOC")
df$tarla <- ifelse(df$tarla == "pre","Pre-Tarla","Tarla")
df$treatment_status <- ifelse(df$treatment_status == "Naive", df$treatment_status, ifelse(df$tarla == "Pre-Tarla" | is.na(df$tarla), df$treatment_status,"Tarla"))

df %>% 
  count(subtype,treatment_status) %>% 
  group_by(treatment_status) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n/total)*100) %>% 
  select(treatment_status,subtype,freq) %>% 
  arrange(treatment_status)
