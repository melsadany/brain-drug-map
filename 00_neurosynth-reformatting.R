################################################################################
#                    rewrite neurosynth data in needed format                  #
################################################################################
rm(list = ls())
gc()
source("/wdata/msmuhammad/msmuhammad-source.R")
#################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/brain-drug-map"
setwd(project.dir)
#################################################################################
ns.dir <- "/sdata/neurosynth/neurosynth-data"
coords <- read_tsv(paste0(ns.dir, "/data-neurosynth_version-7_coordinates.tsv.gz"))
metadata <- read_tsv(paste0(ns.dir, "/data-neurosynth_version-7_metadata.tsv.gz"))
features <- read_csv(paste0(ns.dir, "/features-combined-v7.csv"))
metadata2 <- inner_join(coords, metadata)
write_tsv(metadata2, "data/neurosynth/metadata-v7.tsv")
#################################################################################
# make a dataframe of xyz coordinates as rows and columns to be terms
# make sure to keep xyz coords from MNI space only

# drop features that are not chr
features.meta <- data.frame(f = colnames(features)) %>%
  mutate(sum = colSums(features)) %>%
  mutate(chr = c(rep(F, 60), rep(T, ncol(features)-60))) %>%
  filter(chr == T)
# round the xyz position
all <- left_join(coords%>%select(id,x,y,z),
                 cbind(metadata %>% select(id, ref_space = space), features%>%select(features.meta$f))) %>%
  filter(ref_space == "MNI") %>%
  mutate_at(.vars = vars(x,y,z), .funs = function(x) round(x))
gc()
# aggregate the values to get their mean based on the xyz position
t <- as.data.table(all)
aggregated_data <- t[, lapply(.SD, mean, na.rm = TRUE), 
                              by = .(x, y, z), .SDcols = features.meta$f]
gc()
# save
write_tsv(aggregated_data, "data/neurosynth/features-data-by-xyz-in-MNI-only.tsv")
pdssave(aggregated_data, file = "data/neurosynth/features-data-by-xyz-in-MNI-only.rds")
#################################################################################