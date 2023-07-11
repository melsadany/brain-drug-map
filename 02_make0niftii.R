################################################################################
#                     make a nifti of predictions for a drug                   #
################################################################################
rm(list = ls())
gc()
# source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
library(tidyverse)
library(oro.nifti)
library(neurobase, lib.loc = "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/miniconda3/envs/ENA/lib/R/library")
library(fslr, lib.loc = "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/miniconda3/envs/ENA/lib/R/library")
#################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/brain-drug-map"
setwd(project.dir)

#################################################################################
# use annot for m1, and whole for m2, m3
# mni <- readNIfTI("/Dedicated/jmichaelson-wdata/msmuhammad/refs/mni_icbm152_nlin_sym_09c_CerebrA_nifti/mni_icbm152_CerebrA_tal_nlin_sym_09c.nii", reorient = F)
mni <- readNIfTI("/Dedicated/jmichaelson-wdata/msmuhammad/refs/mni_icbm152_nl_VI_nifti/icbm_avg_152_t1_tal_nlin_symmetric_VI_brain.nii.gz", reorient = F)

atlas.annot <- read_rds("../../refs/mni_icbm152_nlin_sym_09c_CerebrA_nifti/annotated-xyz.rds")
atlas.mni <- read_rds("../../refs/mni_icbm152_nl_VI_nifti/brain-pos.rds")

# origin <- c(97,133,79)
# origin <- c(91,127,73)

m1.predictions <- read_rds("data/model-derivatives/01_m1-glm-lincs-predictions.rds")
m2.predictions <- read_rds("data/model-derivatives/01_m2-glm-lincs-predictions.rds")
m3.predictions <- read_rds("data/model-derivatives/01_m3-glm-lincs-predictions.rds")


glm.predictions <- m2.predictions

meph <- inner_join(glm.predictions %>%
                    select(starts_with("mni"), corr = methylphenidate),
                  # atlas.annot)
                  atlas.mni)

mni2 <- mni
dim(mni)
mni2[as.matrix(meph[,5:7])] <- meph$corr

### Writing out pct without any datatype change - still UINT8
pct <- mni2 / max(mni2)
pct <- cal_img(pct)
pct@bitpix
hist(pct)
hist(mni2)
pct = drop_img_dim(pct)
pct = zero_trans(pct)
pct = onefile(pct)
tfile = "data/01_m2-glm-lincs-mph-predictions"
writeNIfTI(pct, filename = tfile, verbose=TRUE)
pct = datatyper(pct, type_string = "FLOAT32")
writeNIfTI(pct, filename = tfile)
pct3 = readNIfTI(paste0(tfile, ".nii.gz"))
hist(pct3, main="Writing out with changing of datatype/bitpix")
########################################################################################33
# DL

m4.predictions <- read_rds("data/model-derivatives/01_m4-dl-cmap-predictions.rds")
dl.predictions <- m4.predictions

med <- inner_join(dl.predictions %>%
                     select(starts_with("mni"), 
                            # corr = methylphenidate
                            corr = sumatriptan),
                   atlas.mni)
mni2 <- mni
dim(mni)
mni2[as.matrix(med[,5:7])] <- med$corr
tfile = "data/01_m4-dl-cmap-sumatriptan-predictions"
writeNIfTI(mni2, filename = tfile, verbose=TRUE)
