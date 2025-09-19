source("source/sclc_cytof_functions.R")
################################################################################
# Output all CTCs
################################################################################
ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

fcs <- CATALYST::sce2fcs(ctcs, split_by = "collection_id")

write.flowSet(fcs, outdir = "data/output_fcs_files")

################################################################################
# Output Subtypes
################################################################################
fcs <- CATALYST::sce2fcs(ctcs, split_by = "subtype")

write.flowSet(fcs, outdir = "data/output_fcs_files")

################################################################################
# Output Naive and Treated files
################################################################################

# Subset to remove post-tarla cells
treatment_status_ctcs <- ctcs[,is.na(ctcs$tarla) | ctcs$tarla != "post"]

fcs <- CATALYST::sce2fcs(treatment_status_ctcs, split_by = "treatment_status")

write.flowSet(fcs, outdir = "data/output_fcs_files")

################################################################################
# Output Pre- and Post-Tarla Files
################################################################################

# Subset to only tarla cells
tarla_ctcs <- ctcs[,!is.na(ctcs$tarla)]

fcs <- CATALYST::sce2fcs(tarla_ctcs, split_by = "tarla")

write.flowSet(fcs, outdir = "data/output_fcs_files")






