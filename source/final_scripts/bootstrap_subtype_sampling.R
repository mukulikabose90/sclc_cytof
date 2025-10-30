source("source/sclc_cytof_functions.R")

set.seed(42)

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

cluster_colors <- c("#dd4b33", "#F1FAEE", "#A8DADC", "#457B9D")

status_colors <- c("#E63946","#457B9D","#A8DADC")


subtype_status_colors <- c("#b54a3a","#dd4b33","#E6EAE3", "#F1FAEE","#8FB6B8","#A8DADC","#3D617C","#3D617C")

subtype_status_colors <- c(
  "#b54a3a", "#dd4b33", "#f08070",  # shades of red
  "#8FB6B8", "#A8DADC", "#D0EBEB",  # shades of teal
  "#3D617C", "#5A7B95", "#7A99B0",  # shades of blue
  "#E6EAE3", "#F1FAEE", "#FFFFF0"   # shades of off-white/green
)

all_data <- ctcs@colData %>% 
  as.data.frame()

all_data$treatment_status <- ifelse(all_data$treatment_status == "naive","Naive","CTX ± ICI")

all_data$treatment_status <- ifelse(is.na(all_data$tarla), all_data$treatment_status,
                                    ifelse(all_data$tarla == "pre", all_data$treatment_status, "Tarla"))


################################################################################
# Calculate mean subtype proportions
# Naive vs SOC
################################################################################
n_cells <- 30

# data_df <- all_data %>% 
#   filter(treatment_status %in% c("Naive","SOC"))

# Repeat 1000 times
resamples <- map(1:1000, ~ {
  all_data %>% 
    group_by(collection_id) %>%
    filter(n() >= n_cells) %>%
    slice_sample(n = n_cells) %>%
    count(subtype,treatment_status) %>%                   # count subtypes across ALL patients
    mutate(prop = n / sum(n), 
           iteration = .x) %>% 
    ungroup()
})

# Combine
df_resampled <- bind_rows(resamples)

df_resampled$treatment_status <- factor(df_resampled$treatment_status, levels = c("Naive","CTX ± ICI","Tarla"))
df_resampled$subtype <- factor(df_resampled$subtype, levels = c("A","N","P","Mes"))
df_resampled$treatment_subtype <- paste0(df_resampled$treatment_status,"_",df_resampled$subtype)
df_resampled$treatment_subtype <- factor(df_resampled$treatment_subtype, levels=c("Naive_A","SOC_A","Tarla_A", "Naive_N","SOC_N","Tarla_N",
                                                                                  "Naive_P","SOC_P","Tarla_P","Naive_Mes","SOC_Mes","Tarla_Mes"))


stat.test <- df_resampled %>%
  group_by(subtype) %>%
  wilcox_test(prop ~ treatment_status) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance() %>% 
  add_xy_position(x = "treatment_status")

p <- ggviolin(df_resampled, x="treatment_status",y="prop",fill="treatment_status",draw_quantiles = 0.5)+
  facet_wrap(~subtype, scales="free_x",nrow=1)+
  stat_pvalue_manual(stat.test, label = "p.adj.signif",size=5,tip.length = 0)+
  scale_fill_manual(values = status_colors)+
  ylim(0,NA)+
  scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1))+
  labs(y="Proportion",
       x= "Treatment Status")+
  theme(axis.title = element_text(size=20),
        axis.text.x = element_text(size=16,angle=45,hjust=1,vjust=1),
        strip.text = element_text(face = "bold", size=16), 
        strip.background = element_blank())+
  rremove("legend")

p

tiff(glue("figures/resampling_treatment_status_proportion_violinplot.tiff"), width=300,height=200, units = "mm", res=600)
print(p)
dev.off()
