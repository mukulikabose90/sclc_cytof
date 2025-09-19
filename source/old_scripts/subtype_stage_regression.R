

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

ctcs <- ctcs[,!is.na(ctcs$subtype)]

y <- assay(ctcs, "exprs")

df <- data.frame(t(y), colData(ctcs), check.names = FALSE)

temp <- df %>% 
  select(collection_id,stage) %>% 
  dplyr::filter(!is.na(stage)) %>% 
  mutate(y = ifelse(stage == "extended", 1, 0))

cluster_prop_df <- as.data.frame(colData(ctcs)) %>%
  group_by(collection_id,subtype) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))%>% 
  group_by(collection_id) %>% 
  mutate(total = sum(n)) %>% 
  select(collection_id,subtype,freq) %>% 
  pivot_wider(names_from = subtype, values_from = freq)

cluster_prop_df[is.na(cluster_prop_df)] <- 0


data <- merge(temp,cluster_prop_df, by="collection_id",all=T) %>% 
  distinct() %>% 
  dplyr::filter(collection_id != "SC370-3")


model <- lm("y ~ A + I + N + P", data=data)




summary(model)
