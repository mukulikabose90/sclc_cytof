source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################
# Calculate mean subtype proportions
################################################################################
# Repeat 100 times
resamples <- map(1:1000, ~ {
  ctcs@colData %>%
    as.data.frame() %>% 
    group_by(patient_id) %>%
    filter(n() >= 10) %>%
    slice_sample(n = 10) %>%
    count(subtype) %>%                   # count subtypes across ALL patients
    mutate(prop = n / sum(n), 
           iteration = .x) %>% 
    ungroup()
})

# Combine
df_resampled <- bind_rows(resamples)

# Global consensus with variability
global_consensus <- df_resampled %>%
  group_by(subtype) %>%
  summarise(
    mean_prop = mean(prop),
    sd_prop   = sd(prop),
    lower_ci  = quantile(prop, 0.025),
    upper_ci  = quantile(prop, 0.975),
    .groups = "drop"
  )

global_consensus


ggplot(global_consensus, aes(x = subtype, y = mean_prop)) +
  geom_point(size = 3, color = "red") +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
  ylab("Proportion") +
  theme_minimal()

ggplot(global_consensus, aes(x = subtype, y = mean_prop)) +
  geom_point(size = 3, color = "blue") +
  geom_errorbar(aes(ymin = mean_prop - sd_prop, ymax = mean_prop + sd_prop), width = 0.2) +
  ylab("Proportion") +
  theme_minimal()

ggplot(global_consensus, aes(x = subtype, y = mean_prop, fill = subtype)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
  ylab("Proportion") +
  theme_minimal() +
  theme(legend.position = "none")

ggplot(df_resampled, aes(x = subtype, y = prop, fill = subtype)) +
  # geom_violin(trim = FALSE, alpha = 0.4) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  ylab("Proportion per iteration") +
  theme_minimal() +
  theme(legend.position = "none")


################################################################################
# Calculate mean treatment status proportions
################################################################################

df <- ctcs@colData %>%
  as.data.frame()

df$treatment_status <- ifelse(df$treatment_status == "naive","Naive","SOC")
df$tarla <- ifelse(df$tarla == "pre","Pre-Tarla","Tarla")

df$treatment_status <- ifelse(df$treatment_status == "Naive", df$treatment_status, ifelse(df$tarla == "Pre-Tarla" | is.na(df$tarla), df$treatment_status,"Tarla"))



# Repeat 100 times
resamples <- map(1:100, ~ {
  df %>% 
    group_by(patient_id) %>%
    filter(n() >= 10) %>%
    slice_sample(n = 10) %>%
    count(treatment_status) %>%                   # count subtypes across ALL patients
    mutate(prop = n / sum(n), 
           iteration = .x) %>% 
    ungroup()
})

# Combine
df_resampled <- bind_rows(resamples)

# Global consensus with variability
global_consensus <- df_resampled %>%
  group_by(treatment_status) %>%
  summarise(
    mean_prop = mean(prop),
    sd_prop   = sd(prop),
    lower_ci  = quantile(prop, 0.025),
    upper_ci  = quantile(prop, 0.975),
    .groups = "drop"
  )

global_consensus



ggplot(global_consensus, aes(x = treatment_status, y = mean_prop)) +
  geom_point(size = 3, color = "red") +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
  ylab("Proportion") +
  theme_minimal()
