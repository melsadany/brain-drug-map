################################################################################
#               predicting activity of drugs in brain xyz positions            #
################################################################################
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
# library(tidyverse)
library(oro.nifti, lib.loc = "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/miniconda3/envs/ENA/lib/R/library")
library(tensorflow)
use_condaenv("EDL")
library(keras)
k_clear_session()
set.seed(123)
# pload <- function(fname,envir=.GlobalEnv){
#   con <- pipe(paste("/Dedicated/jmichaelson-wdata/msmuhammad/workbench/pixz -d <",fname),"rb")
#   load(con,envir=envir); close(con)
# }
#################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/brain-drug-map"
setwd(project.dir)
#################################################################################
load("data/exp/exp-corr-sex-age-ethn-donor-PMI_w-annot.rda")
rm(avg.brain.exp.ls)
gc()
allen.annotated <- cbind(annot.ls, corr.brain.exp.ls)
corr.brain.exp <- allen.annotated%>%select(colnames(corr.brain.exp.ls))
annot.ls.filt <- allen.annotated%>%select(starts_with("mni")&ends_with("allen"),region)
#################################################################################
#################################################################################
pload("/Dedicated/jmichaelson-sdata/CMAP/LINCS/small_mol_metadata_w_BBB_FDA.Rdata.pxz")
bbb.drugs <- smol_meta %>%
  filter(FDA == T, BBB_perm>0.5)
drug.sig <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/LINCS/cmap.of.int_v2.rds") %>%
  as.data.frame() %>%
  select(any_of(colnames(corr.brain.exp.ls))) %>%
  t() %>%
  as.data.frame() %>%
  select(any_of(bbb.drugs$pert_name))
# drug.sig <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/LINCS/LINCS-combined-corrected-14.rds") %>%
#   arrange(symbol) %>%
#   column_to_rownames("symbol")
# drug.sig <- drug.sig[colnames(corr.brain.exp),]
corr.brain.exp <- corr.brain.exp[,rownames(drug.sig)]
table(rownames(drug.sig) == colnames(corr.brain.exp))
#################################################################################
# compute correlations with drugs
drug.corr <- cor(t(corr.brain.exp), drug.sig, method = "spearman")
drug.corr.annot <- cbind(annot.ls.filt, drug.corr)
# transformed the correlatios to be from -1 to 1 per drug
trans.drug.corr <- apply(drug.corr, MARGIN = 2, FUN = function(x) ((x-min(x))/(max(x)-min(x)))* 2 - 1)
#################################################################################
all <- cbind(annot.ls.filt, trans.drug.corr)
# ohe.all <-model.matrix(~ region - 1, 
#                        all%>%mutate(region = as.factor(region)))
# colnames(ohe.all) <- sub("region", "", colnames(ohe.all))
# all.mod <- cbind(all, ohe.all) 
all.mod <- all
x.train <- all.mod %>%
  select(c(starts_with("mni")
           ))
y.train <- all.mod %>%
  select(colnames(drug.sig))
# hist(x.train$mni_x_allen)
# hist(x.train$mni_y_allen)
# hist(x.train$mni_z_allen)
# hist(y.train$methylphenidate)
#################################################################################
# Define the neural network architecture
rm(model)
k_clear_session()

model <- keras_model_sequential(input_shape = c(dim(x.train)[2])) %>%
  layer_dense(units = 128, activation = "relu",
              # use_bias = TRUE,
              # kernel_regularizer = regularizer_l2(1e-6),
              ) %>%
  layer_dense(units = 64, activation = "relu",
              # use_bias = TRUE,
              # kernel_regularizer = regularizer_l2(1e-6),
  ) %>%
  layer_dense(units = 64, activation = "relu",
              # use_bias = TRUE,
              # kernel_regularizer = regularizer_l2(1e-6),
  ) %>%
  # layer_dropout(0.2) %>%
  layer_dense(units = 32, activation = "relu",
              use_bias = TRUE,
              kernel_regularizer = regularizer_l2(1e-6),
              ) %>%
  layer_dense(units = ncol(y.train), activation = "tanh", 
              kernel_regularizer=regularizer_l2(1e-6),
              use_bias = T)
#################################################################################
#################################################################################
# Compile the model
model %>% compile(
  loss = "mse",
  optimizer = optimizer_rmsprop(learning_rate = 0.01),
  metrics = c("mse", "mae")
)
# summary(model)
# Train the model
history <- model %>% fit(
  x = x.train%>%as.matrix(), 
  y = y.train%>%as.matrix(),
  epochs = 10,
  batch_size = 64, 
  validation_split = 0.01
)
save_model_tf(model, "data/model-derivatives/model-cmap-FDA-BBB0.5-071123")
# mni <- readNIfTI("/Dedicated/jmichaelson-wdata/msmuhammad/refs/mni_icbm152_nlin_sym_09c_CerebrA_nifti/mni_icbm152_CerebrA_tal_nlin_sym_09c.nii", reorient = F)
# origin <- c(97,133,79)
# mni.hit <- which(mni>0, arr.ind = T) %>%
#   as.data.frame() %>%
#   mutate(mni_x = dim1-origin[1]) %>%
#   mutate(mni_y = dim2-origin[2]) %>%
#   mutate(mni_z = dim3-origin[3])
# mni.hit.labeles <- mni.hit
mni <- readNIfTI("/Dedicated/jmichaelson-wdata/msmuhammad/refs/fsl-data/standard/MNI152_T1_1mm_brain.nii.gz", reorient = F)
origin <- c(91,127,73)
mni.hit <- which(mni>0, arr.ind = T) %>%
  as.data.frame() %>%
  mutate(mni_x = -1*(dim1-origin[1])) %>%
  mutate(mni_y = dim2-origin[2]) %>%
  mutate(mni_z = dim3-origin[3])
mni.hit.whole <- mni.hit

# sample <- sample(1:nrow(mni.hit), size = 200000)
# mni.hit <- mni.hit[sample,]
# mni.hit.whole <- mni.hit.whole[sample,]

predictions <- predict(model, mni.hit%>%select(starts_with("mni"))%>%as.matrix())
colnames(predictions) <- colnames(y.train)
predictions <- cbind(mni.hit.whole, predictions)
write_rds(predictions, "data/model-derivatives/dl-predictions-cmap_whole-brain_FDA-BBB0.5_071123.rds", compress="gz")
##small subset
small <- predictions %>% select(1:6, methylphenidate, venlafaxine, guanfacine, sertraline, clozapine, bupropion, aspirin, atomoxetine, caffeine, citalopram, clonidine, clonazepam, diazepam, fluoxetine, fluvoxamine, haloperidol, histamine, ibuprofen, imatinib, imipramine, ketamine, ketoprofen, levetiracetam, loxapine, minoxidil, naproxen, nicotine, paroxetine, risperidone, tramadol, trazodone, venlafaxine, verapamil, zolpidem)
write_rds(small, "data/model-derivatives/dl-predictions-cmap-subset_whole-brain_FDA-BBB0.5_071123.rds", compress="gz")
#################################################################################
# make nifti per drug
drugs <- colnames(drug.sig)
library(doMC);registerDoMC(6)
foreach::foreach(i=1:length(drugs)) %dopar% {
  med <- predictions %>%
    select(1:6, i+6) %>%
    rename(corr = 7)
  mni2 <- mni
  mni2[as.matrix(med[,1:3])] <- med$corr
  ### Writing out pct without any datatype change - still UINT8
  pct <- mni2 / max(mni2)
  pct <- cal_img(pct)
  # pct@bitpix
  pct = drop_img_dim(pct)
  pct = zero_trans(pct)
  pct = onefile(pct)
  tfile = paste0("data/maps/dl-cmap_whole-brain_FDA-BBB0.5_071123_", drugs[i])
  writeNIfTI(pct, filename = tfile, verbose=TRUE)
  pct = datatyper(pct, type_string = "FLOAT32")
  writeNIfTI(pct, filename = tfile)
}
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
# try another approach
# predict gene expression vector for all brain regions instead of drug corr
rm(list = ls())
gc()
library(tidyverse, lib.loc = "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/miniconda3/envs/EDL/lib/R/library")
library(oro.nifti, lib.loc = "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/miniconda3/envs/ENA/lib/R/library")
library(tensorflow)
library(doMC)
use_condaenv("EDL")
library(keras)
k_clear_session()
set.seed(123)
#################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/brain-drug-map"
setwd(project.dir)
#################################################################################
load("data/exp/exp-corr-sex-age-ethn-donor-PMI_w-annot.rda")
# cmap.of.int <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/LINCS/cmap.of.int.rds")
# genes <- colnames(cmap.of.int)
#################################################################################
xyz <- annot.ls[,paste0("mni_", c("x", "y", "z"))] %>% as.matrix()
x.train <- xyz
y.train <- corr.brain.exp.ls %>% as.matrix()
dim(y.train)
####
# get stats about genes
df <- cbind(apply(y.train, MARGIN = 2, function(x) mean(x)) %>% as.data.frame(),
            apply(y.train, MARGIN = 2, function(x) var(x)) %>% as.data.frame()) %>% 
  rownames_to_column("gene") %>%
  rename(mean = 2, var = 3)
df %>% 
  ggplot(aes(x=mean, y=var)) +
  geom_point(size=0.3) +
  geom_vline(xintercept = 0)+geom_smooth()
# use gene expression to predict region, and if it's significant, keep it
registerDoMC(cores = 6)
r2 = foreach(i=1:ncol(y.train), .combine = rbind) %dopar% summary(lm(y.train[,i]~annot.ls$structure_acronym))$adj.r.squared
rownames(r2) = colnames(y.train)
# r2 = unlist(r2)
write_rds(r2, "data/exp/r2-of-corrected-gene-predicting-structure.rds", compress = "gz")
r2 <- read_rds("data/exp/r2-of-corrected-gene-predicting-structure.rds")
# filter y.train to keep genes of interest
# y.train <- y.train[,df$gene[df$var>0.6]]
y.train <- y.train[,rownames(r2)[r2[,1]>0.5]]
gc()
#################################################################################
# model
rm(model)
gc()
k_clear_session()
model <- keras_model_sequential(input_shape = c(dim(xyz)[2])) %>%
  # layer_dense(units=512,activation="tanh",
  #             use_bias=T,
  # ) %>%
  # # layer_dropout(0.2) %>%
  # # layer_dense(units=512,activation="relu",
  # #             use_bias=T,
  # # ) %>%
  # layer_dense(units=256,activation="tanh",
  #             use_bias=T,
  # ) %>%
  # # layer_dense(units=64,activation="relu",
  # #             use_bias=T,
  # # ) %>%
  # # layer_dropout(0.33) %>%
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
  layer_dense(units=256,activation="linear",
              kernel_initializer=initializer_orthogonal(),
              use_bias=F,
  ) %>%
  layer_dropout(0.33) %>%
  # layer_dense(units=ncol(y.train),activation="tanh",
  #             kernel_regularizer=regularizer_l1(1e-4),
  #             #	kernel_constraint=constraint_nonneg(),
  #             use_bias=T)
  layer_dense(units=ncol(y.train),activation="linear",
              kernel_regularizer=regularizer_l1(1e-4),
              use_bias=T)
#################################################################################
# Compile the model
model %>% compile(
  loss = "mse",
  optimizer = optimizer_adam(5e-4,clipnorm=0.1),
  weighted_metrics = 'mae',
  loss_weights=list(c(1))
)
# summary(model)
# Train the model
history <- model %>% fit(
  x = x.train,
  y = y.train,
  epochs = 1000,
  batch_size = 512,
  validation_split = 0.01,
  verbose=1,
  callbacks = list(callback_tensorboard())
  # validation_data = list(x.val%>%as.matrix(), y.val%>%as.matrix())
)
save_model_tf(model, "data/model-derivatives/exp-prediction-model-r205-081723")
mni <- readNIfTI("/Dedicated/jmichaelson-wdata/msmuhammad/refs/fsl-data/standard/MNI152_T1_1mm_brain.nii.gz", reorient = F)
origin <- c(91,127,73)
mni.hit <- which(mni>0, arr.ind = T) %>%
  as.data.frame() %>%
  mutate(mni_x = -1*(dim1-origin[1])) %>%
  mutate(mni_y = dim2-origin[2]) %>%
  mutate(mni_z = dim3-origin[3])
mni.hit.whole <- mni.hit
predictions <- predict(model, mni.hit.whole[,4:6])
colnames(predictions) <- rownames(r2)[r2[,1]>0.5]
pdssave(cbind(mni.pos, predictions), file = "data/model-derivatives/predicted-exp-whole-brain-081723.rds")
#################################################################################
# compute correlation with drug signatures
predictions <- pdsload("data/model-derivatives/predicted-exp-whole-brain-081723.rds.pxz")[,-c("x", "y", "z")]
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
#################################################################################
# compute correlations with drugs
drug.corr <- cor(t(predictions), drug.sig, method = "spearman")
# drug.annot <- cbind(mni.pos, drug.corr)
pdssave(cbind(mni.pos, drug.corr), file = "data/model-derivatives/all-drug-map-predicted-exp-whole-brain-081723.rds")
# pdssave(drug.annot[,c("x", "y", "z", "methylphenidate", "")], "data/model-derivatives/predicted-all-drug-map-exp-whole-brain-081723.rds")
#################################################################################
#################################################################################
small <- cbind(mni.pos, drug.corr) %>% as.data.frame()%>% select(x,y,z, 
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
