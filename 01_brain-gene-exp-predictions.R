################################################################################
#                predicting gene expression in brain xyz positions             #
################################################################################
rm(list = ls())
gc()
.libPaths("/old_Users/msmuhammad/workbench/miniconda3/envs/EDL/lib/R/library")
library(tidyverse, lib.loc = "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/miniconda3/envs/EDL/lib/R/library")
library(oro.nifti, lib.loc = "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/miniconda3/envs/ENA/lib/R/library")
library(tensorflow)
library(doMC)
use_condaenv("EDL")
library(keras)
k_clear_session()
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
load("data/exp/exp-corr-sex-age-ethn-donor-PMI_w-annot.rda")
#################################################################################
xyz <- annot.ls[,paste0("mni_", c("x", "y", "z"))] %>% as.matrix()
x.train <- xyz
y.train <- corr.brain.exp.ls %>% as.matrix()
dim(y.train)
####
r2 <- read_rds("data/exp/r2-of-corrected-gene-predicting-structure.rds")
# filter y.train to keep genes of interest
y.train <- y.train[,rownames(r2)[r2[,1]>0.5]]
gc()
#################################################################################
# model
rm(model)
gc()
k_clear_session()
ipt <- layer_input(shape=c(dim(xyz)[2]))
m1 <- ipt %>%
  layer_dense(units=32,activation="relu",
              use_bias=T,
  ) %>%
  layer_dropout(0.2) %>%
  layer_dense(units=32,activation="relu",
              use_bias=T,
  ) %>%
  layer_dense(units=32,activation="relu",
              use_bias=T,
  ) %>%
  layer_batch_normalization() %>%
  layer_dense(units=128,activation="linear",
              kernel_constraint=constraint_minmaxnorm(0.1,10,rate=0.9,axis=1),
              kernel_regularizer=regularizer_l1(5e-5),
              use_bias=F,
  ) %>%
  layer_dropout(0.33)
out_expr <- m1 %>%
  layer_dense(units=ncol(y.train),activation="sigmoid",
              kernel_regularizer=regularizer_l1(1e-5),
              use_bias=T)
out_int <- m1 %>%
  layer_dense(units=1,activation="linear",
              kernel_regularizer=regularizer_l1(1e-5),
              use_bias=T)

model <- keras_model(ipt,list(out_expr,out_int))
#################################################################################
# Compile the model
model %>% compile(
  optimizer = optimizer_adam(1e-3),
  loss = list("mse",'mse'),
  weighted_metrics = list('mae',"mae"),
  loss_weights=list(c(1),c(10))
)
# summary(model)
# Train the model
history <- model %>% fit(
  x = x.train,
  y = y.train,
  epochs = 100,
  batch_size = ncol(y.train),
  validation_split = 0.01,
  verbose=1,
  callbacks = list(callback_tensorboard()),
  sample_weight=list(w0,w0)
  )
w <- get_weights(model)
save_model_tf(model, "data/model-derivatives/gene-exp-prediction-model-r205-092723")

mni <- readNIfTI("/Dedicated/jmichaelson-wdata/msmuhammad/refs/fsl-data/standard/MNI152_T1_1mm_brain.nii.gz", reorient = F)
origin <- c(91,127,73)
mni.hit <- which(mni>0, arr.ind = T) %>%
  as.data.frame() %>%
  mutate(mni_x = -1*(dim1-origin[1])) %>%
  mutate(mni_y = dim2-origin[2]) %>%
  mutate(mni_z = dim3-origin[3])
mni.hit.whole <- mni.hit %>% as.matrix()
predictions <- predict(model, mni.hit.whole[,4:6])
colnames(predictions) <- rownames(r2)[r2[,1]>0.5]
pdssave(cbind(mni.hit.whole, predictions), file = "data/model-derivatives/gene-exp-whole-brain-092723.rds")
#################################################################################
