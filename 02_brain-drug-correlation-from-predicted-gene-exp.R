################################################################################
#                 compute drug correlation in brain xyz positions              #
################################################################################
rm(list = ls())
gc()
library(tidyverse, lib.loc = "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/miniconda3/envs/EDL/lib/R/library")
library(oro.nifti, lib.loc = "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/miniconda3/envs/ENA/lib/R/library")
library(doMC)
set.seed(123)
pdsload <- function(fname,envir=.GlobalEnv){
  con <- pipe(paste("/Dedicated/jmichaelson-wdata/msmuhammad/workbench/pixz -d <",fname),"rb")
  return(readRDS(con)); close(con)
}
pdssave <- function(...,file){  
  con = pipe(paste("/Dedicated/jmichaelson-wdata/msmuhammad/workbench/pixz -2 -q 80 -f 3 > ",file,".pxz",sep=""),"wb") 
  saveRDS(...,file=con); close(con) 
} 
pload <- function(fname,envir=.GlobalEnv){
  con <- pipe(paste("/Dedicated/jmichaelson-wdata/msmuhammad/workbench/pixz -d <",fname),"rb")
  load(con,envir=envir); close(con)
}
#################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/brain-drug-map"
setwd(project.dir)
#################################################################################
# in case of running out of memory after saving the predictions file, only do the correlation again
#################################################################################
# compute correlation with drug signatures
# predictions <- pdsload("data/model-derivatives/gene-exp-whole-brain-082223.rds.pxz")
predictions <- pdsload("data/model-derivatives/gene-exp-whole-brain-drug-genes-082223.rds.pxz")
pload("/Dedicated/jmichaelson-sdata/CMAP/LINCS/small_mol_metadata_w_BBB_FDA.Rdata.pxz")
bbb.drugs <- smol_meta %>%
  filter(FDA == T, BBB_perm>0.5)
drug.sig <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/LINCS/cmap.of.int_v2.rds") %>%
  as.data.frame() %>%
  select(any_of(colnames(predictions))) %>%
  t() %>%
  as.data.frame() %>%
  select(any_of(bbb.drugs$pert_name))
predictions <- scale(predictions[,rownames(drug.sig)], scale = T, center = T)
table(rownames(drug.sig) == colnames(predictions))
# pdssave(predictions, file = "data/model-derivatives/gene-exp-whole-brain-drug-genes-082223.rds")
#################################################################################
mni <- readNIfTI("/Dedicated/jmichaelson-wdata/msmuhammad/refs/fsl-data/standard/MNI152_T1_1mm_brain.nii.gz", reorient = F)
origin <- c(91,127,73)
mni.hit <- which(mni>0, arr.ind = T) %>%
  as.data.frame() %>%
  mutate(mni_x = -1*(dim1-origin[1])) %>%
  mutate(mni_y = dim2-origin[2]) %>%
  mutate(mni_z = dim3-origin[3])
mni.hit.whole <- mni.hit
# compute correlations with drugs
drug.corr <- cor(t(predictions), drug.sig)
pdssave(cbind(mni.hit.whole, drug.corr), file = "data/all-drug-map-predicted-exp-whole-brain-082223.rds")
#################################################################################
small <- cbind(mni.hit.whole, drug.corr) %>% as.data.frame() %>% select(colnames(mni.hit.whole), 
                                                                  methylphenidate, venlafaxine, guanfacine, 
                                                                  sertraline, clozapine, bupropion, aspirin, 
                                                                  atomoxetine, caffeine, citalopram, clonidine, 
                                                                  clonazepam, diazepam, fluoxetine, fluvoxamine, 
                                                                  haloperidol, histamine, ibuprofen, imatinib, 
                                                                  imipramine, ketamine, ketoprofen, levetiracetam, 
                                                                  loxapine, minoxidil, naproxen, nicotine, paroxetine, 
                                                                  risperidone, tramadol, trazodone, venlafaxine, 
                                                                  verapamil, zolpidem)
pdssave(small, file = "data/subsetof-drug-map-predicted-exp-whole-brain-082223.rds")
#################################################################################
#################################################################################
