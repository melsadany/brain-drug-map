


rm(list = ls())
gc()
library(tidyverse, lib.loc = "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/miniconda3/envs/EDL/lib/R/library")
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
#################################################################################
# compute correlation with drug signatures
# predictions <- pdsload("data/model-derivatives/predicted-exp-whole-brain-081723.rds.pxz")
predictions <- pdsload("data/predicted-exp-whole-brain-drug-genes-081723.rds.pxz")
pload("/Dedicated/jmichaelson-sdata/CMAP/LINCS/small_mol_metadata_w_BBB_FDA.Rdata.pxz")
bbb.drugs <- smol_meta %>%
  filter(FDA == T, BBB_perm>0.5)
drug.sig <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/LINCS/cmap.of.int_v2.rds") %>%
  as.data.frame() %>%
  select(any_of(colnames(predictions))) %>%
  t() %>%
  as.data.frame() %>%
  select(any_of(bbb.drugs$pert_name))
predictions <- predictions[,rownames(drug.sig)]
table(rownames(drug.sig) == colnames(predictions))
# pdssave(predictions, file = "data/predicted-exp-whole-brain-drug-genes-081723.rds")
#################################################################################
library(oro.nifti, lib.loc = "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/miniconda3/envs/ENA/lib/R/library")
mni <- readNIfTI("/Dedicated/jmichaelson-wdata/msmuhammad/refs/mni_icbm152_nlin_sym_09c_CerebrA_nifti/mni_icbm152_CerebrA_tal_nlin_sym_09c.nii", reorient = F)
mni.pos <- which(mni>0, arr.ind = T)
origin <- c(97,133,79)
mni.pos <- mni.pos %>%
  as.data.frame() %>%
  mutate(x=dim1-origin[1]) %>%
  mutate(y=dim2-origin[2]) %>%
  mutate(z=dim3-origin[3]) %>%
  select(x,y,z) %>% as.matrix()

# compute correlations with drugs
drug.corr <- cor(t(predictions), drug.sig, method = "spearman")
# drug.annot <- cbind(mni.pos, drug.corr)
pdssave(cbind(mni.pos, drug.corr), file = "data/all-drug-map-predicted-exp-whole-brain-081723.rds")
# pdssave(drug.annot[,c("x", "y", "z", "methylphenidate", "")], "data/model-derivatives/predicted-all-drug-map-exp-whole-brain-081723.rds")
#################################################################################
small <- cbind(mni.pos, drug.corr) %>% as.data.frame() %>% select(x,y,z, 
                                              methylphenidate, venlafaxine, guanfacine, 
                                              sertraline, clozapine, bupropion, aspirin, 
                                              atomoxetine, caffeine, citalopram, clonidine, 
                                              clonazepam, diazepam, fluoxetine, fluvoxamine, 
                                              haloperidol, histamine, ibuprofen, imatinib, 
                                              imipramine, ketamine, ketoprofen, levetiracetam, 
                                              loxapine, minoxidil, naproxen, nicotine, paroxetine, 
                                              risperidone, tramadol, trazodone, venlafaxine, 
                                              verapamil, zolpidem)
pdssave(small, file = "data/subsetof-drug-map-predicted-exp-whole-brain-081723.rds")
#################################################################################
#################################################################################
#################################################################################
#################################################################################
# independent
library(oro.nifti, lib.loc = "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/miniconda3/envs/ENA/lib/R/library")

mni <- readNIfTI("/Dedicated/jmichaelson-wdata/msmuhammad/refs/mni_icbm152_nlin_sym_09c_CerebrA_nifti/mni_icbm152_CerebrA_tal_nlin_sym_09c.nii", reorient = F)
mni.pos <- which(mni>0, arr.ind = T)
origin <- c(97,133,79)
mni.pos <- mni.pos %>%
  as.data.frame() %>%
  mutate(x=dim1-origin[1]) %>%
  mutate(y=dim2-origin[2]) %>%
  mutate(z=dim3-origin[3])

pload("/Dedicated/jmichaelson-sdata/CMAP/LINCS/small_mol_metadata_w_BBB_FDA.Rdata.pxz")
bbb.drugs <- smol_meta %>%
  filter(FDA == T, BBB_perm>0.5)
drug.sig <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/LINCS/cmap.of.int_v2.rds") %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  select(any_of(bbb.drugs$pert_name))
drug.corr <- pdsload(fname = "data/all-drug-map-predicted-exp-whole-brain-081723.rds.pxz") %>% as.data.frame()
# make nifti per drug
drugs <- colnames(drug.sig)
library(doMC);registerDoMC(6)
foreach::foreach(i=1:length(drugs)) %dopar% {
  d <- colnames(drug.sig)[i]
  med <- inner_join(mni.pos,
                    drug.corr %>%
                      select(x,y,z, d) %>%
                      rename(corr = 4))
  mni2 <- mni
  mni2[as.matrix(med[,1:3])] <- med$corr
  ### Writing out pct without any datatype change - still UINT8
  pct <- mni2 / max(mni2)
  pct <- cal_img(pct)
  # pct@bitpix
  pct = drop_img_dim(pct)
  pct = zero_trans(pct)
  pct = onefile(pct)
  tfile = paste0("data/maps/model-081723/all/", drugs[i])
  writeNIfTI(pct, filename = tfile, verbose=TRUE)
  pct = datatyper(pct, type_string = "FLOAT32")
  writeNIfTI(pct, filename = tfile)
}
#################################################################################
