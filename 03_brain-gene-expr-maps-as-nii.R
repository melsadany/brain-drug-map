################################################################################
#                  save drug correlations as brain nifti images                #
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
library(neurobase, lib.loc = "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/miniconda3/envs/ENA/lib/R/library")
library(fslr, lib.loc = "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/miniconda3/envs/ENA/lib/R/library")
#################################################################################
#################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/brain-drug-map"
setwd(project.dir)
#################################################################################
#################################################################################
mni <- readNIfTI("/Dedicated/jmichaelson-wdata/msmuhammad/refs/fsl-data/standard/MNI152_T1_1mm_brain.nii.gz", reorient = F)
origin <- c(91,127,73)
mni.hit <- which(mni>0, arr.ind = T) %>%
  as.data.frame() %>%
  mutate(mni_x = -1*(dim1-origin[1])) %>%
  mutate(mni_y = dim2-origin[2]) %>%
  mutate(mni_z = dim3-origin[3])
mni.hit.whole <- mni.hit
#################################################################################
gene.map <- pdsload(fname = "data/model-derivatives/gene-exp-whole-brain-082223.rds.pxz") %>% as.data.frame()
print("done reading the mar")
#################################################################################
# make nifti per gene
genes <- colnames(gene.map)
library(doMC);registerDoMC(6)
system("mkdir -p data/maps/model-082223/gene-exp")
foreach::foreach(i=1:length(genes)) %dopar% {
  d <- colnames(gene.map)[i]
  med <- inner_join(mni.hit.whole,
                    gene.map %>%
                      select(mni_x,mni_y,mni_z, d) %>%
                      rename(corr = 4)) %>%
    mutate(corr = scale(corr, scale = T, center = T))
  mni2 <- mni
  mni2[as.matrix(med[,1:3])] <- med$corr
  ### Writing out pct without any datatype change - still UINT8
  pct <- mni2 / max(mni2)
  pct <- cal_img(pct)
  # pct@bitpix
  pct = drop_img_dim(pct)
  pct = zero_trans(pct)
  pct = onefile(pct)
  tfile = paste0("data/maps/model-082223/gene-exp/", genes[i])
  writeNIfTI(pct, filename = tfile, verbose=TRUE)
  pct = datatyper(pct, type_string = "FLOAT32")
  writeNIfTI(pct, filename = tfile)
  print(paste0("done with gene:", d, ", which is number: ", i, "/", length(genes)))
}


#################################################################################

#################################################################################

#################################################################################