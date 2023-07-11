################################################################################
#              predeicting activity of drugs in brain xyz positions            #
################################################################################
rm(list = ls())
gc()
# source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
library(tidyverse)
library(oro.nifti, lib.loc = "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/miniconda3/envs/ENA/lib/R/library")

library(tensorflow)
use_condaenv("EDL")
library(keras)
k_clear_session()
set.seed(123)
#################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/brain-drug-map"
setwd(project.dir)

#################################################################################
avg.brain.exp.ls <- read_rds("data/exp/combined-aggr-exp-of-all-subjects.rds") %>% as.data.frame()
annot.ls <- read_rds("data/exp/combined-annot-of-all-subjects.rds") %>% as.data.frame()

#################################################################################
# add the atlas labels instead of the allen annotation
atlas.labels <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/refs/mni_icbm152_nlin_sym_09c_CerebrA_nifti/annotated-xyz.rds")
tmp <- left_join(annot.ls %>%
                   rename(mni_x_allen=mni_x)  %>%
                   rename(mni_y_allen=mni_y)  %>%
                   rename(mni_z_allen=mni_z)  %>%
                   mutate(mni_x = round(mni_x_allen))  %>%
                   mutate(mni_y = round(mni_y_allen))  %>%
                   mutate(mni_z = round(mni_z_allen)),
                 atlas.labels)
# write_tsv(tmp, "data/exp/combined-annot-of-all-subjects_w-atlas annot.tsv")
tmp.2 <- cbind(tmp, avg.brain.exp.ls) %>% drop_na(region)
avg.brain.exp.ls.filt <- tmp.2[,29:ncol(tmp.2)]
annot.ls.filt <- tmp.2[,1:28]
#################################################################################
scaled.brain.exp <- scale(avg.brain.exp.ls.filt, center = T, scale = T)
# check if there's a correlation between measured gene exp and sex, age, ethnicity, PMI
summary(lm(data = data.frame(exp = scaled.brain.exp[,100], annot.ls.filt %>% select(donor, sex, age, ethnicity, PMI)),
           exp ~ sex)) #sig for some genes
summary(lm(data = data.frame(exp = scaled.brain.exp[,1], annot.ls.filt %>% select(donor, sex, age, ethnicity, PMI)),
           exp ~ age)) #sig
summary(lm(data = data.frame(exp = scaled.brain.exp[,1], annot.ls.filt %>% select(donor, sex, age, ethnicity, PMI)),
           exp ~ PMI)) #sig
summary(lm(data = data.frame(exp = scaled.brain.exp[,1], annot.ls.filt %>% select(donor, sex, age, ethnicity, PMI)),
           exp ~ donor)) #sig
summary(lm(data = data.frame(exp = scaled.brain.exp[,1], annot.ls.filt %>% select(donor, sex, age, ethnicity, PMI)),
           exp ~ ethnicity)) #sig
summary(lm(data = data.frame(exp = scaled.brain.exp[,1], annot.ls.filt %>% select(donor, sex, age, ethnicity, PMI)),
           exp ~ sex+age+PMI+ethnicity+donor))
# correct avg.exp for each gene based on these covariates
library(doMC)
registerDoMC(cores = 8)
corr.brain.exp.ls <- foreach(i = 1:ncol(avg.brain.exp.ls.filt), .combine = cbind) %dopar% {
  # i=12
  gene <- colnames(avg.brain.exp.ls.filt)[i]
  df <- data.frame(exp = avg.brain.exp.ls.filt[,i], 
                   annot.ls.filt) %>%
    mutate(sc_exp = scale(exp, scale = T, center = T)[,1])
  df <- df%>%
    mutate(corr.exp = residuals(lm(data = df, sc_exp ~ sex+age+PMI+ethnicity)))
  df.2 <- left_join(annot.ls.filt, df) %>%
    select(corr.exp)
  colnames(df.2) <- gene
  return(df.2)
}
save(corr.brain.exp.ls, avg.brain.exp.ls.filt, annot.ls.filt, file = "data/exp/filt-ecp-corr-sex-age-ethn-PMI_w-annot.rda")
# load("data/exp/filt-ecp-corr-sex-age-ethn-PMI_w-annot.rda")
# write_rds(corr.brain.exp.ls, "data/exp/combined-aggr-exp-of-all-subjects_corr-sex-age-donor-ethn-PMI.rds")
# corr.brain.exp.ls <- read_rds("data/exp/combined-aggr-exp-of-all-subjects_corr-sex-age-donor-ethn-PMI.rds")
#################################################################################
# drug.sig <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/LINCS/cmap.of.int.rds") %>% 
#   as.data.frame() %>%
#   select(any_of(colnames(corr.brain.exp.ls))) %>%
#   t()
drug.sig <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/LINCS/LINCS-combined-corrected-14.rds") %>%
  dplyr::filter(symbol %in% colnames(corr.brain.exp.ls)) %>%
  column_to_rownames("symbol")

table(rownames(drug.sig) %in% colnames(corr.brain.exp.ls))

scaled.brain.exp.2 <- scale(corr.brain.exp.ls,scale = T, center = T)
scaled.brain.exp.2 <- scaled.brain.exp.2[,rownames(drug.sig)] %>% t()
all(rownames(scaled.brain.exp.2) == rownames(drug.sig))
scaled.drug.sig <- t(scale(t(drug.sig), scale = T, center = T))

dim(scaled.brain.exp.2)
dim(scaled.drug.sig)
#################################################################################
# compute correlations with drugs
# I transform the drug correlations
drug.corr <- cor(scaled.brain.exp.2, scaled.drug.sig)
drug.corr.annot <- cbind(annot.ls.filt, drug.corr)

# you don't need scaling?
# scaled.drug.corr <- scale(drug.corr, scale = T, center = T)
trans.drug.corr <- drug.corr
for (i in 1:ncol(drug.corr)) {
  d.cor <- drug.corr[,i]
  min.scaled <- min(d.cor)
  max.scaled <- max(d.cor)
  trans.drug.corr[,i] <- (d.cor - min.scaled) / (max.scaled - min.scaled) * 2 - 1
}
#################################################################################
xyz  <- annot.ls.filt
Xs_pre <- xyz %>% 
  select(mni_x_allen, mni_y_allen, mni_z_allen, region)
all <- cbind(Xs_pre, 
             trans.drug.corr
             # drug.corr
             )
ohe.all <-model.matrix(~ region - 1, 
                       all%>%mutate(region = as.factor(region)))
colnames(ohe.all) <- sub("region", "", colnames(ohe.all))
all.mod <- cbind(all, ohe.all) 


# library(caTools)
# index1 <- sample.split(all.mod$region, SplitRatio = 0.8)
# all.2 <- all.mod[index1,]
x.train <- all.mod %>%
# x.train <- all.mod[index1,] %>%
  select(c(starts_with("mni")
           # , colnames(ohe.all)
           ))
y.train <- all.mod %>%
# y.train <- all.mod[index1,] %>%
  select(colnames(drug.sig))
  # select(methylphenidate)

# x.test <- all.mod[!index1,] %>%
#   select(c(starts_with("mni")
#            , colnames(ohe.all)
#            ))
# y.test <- all.mod[!index1,] %>%
#   select(colnames(drug.sig))
# 
hist(x.train$mni_x_allen)
hist(x.train$mni_y_allen)
hist(x.train$mni_z_allen)
hist(y.train$methylphenidate)
#################################################################################
# Define the neural network architecture
rm(model)
k_clear_session()

model <- keras_model_sequential(input_shape = c(dim(x.train)[2])) %>%
  layer_dense(units = 64, activation = "relu",
              # use_bias = TRUE,
              # kernel_regularizer = regularizer_l2(1e-6),
              ) %>%
  # layer_dropout(0.2) %>%
  layer_dense(units = 64, activation = "relu",
              use_bias = TRUE,
              kernel_regularizer = regularizer_l2(1e-6),
              ) %>%
  layer_dense(units = ncol(y.train), activation = "tanh", 
              kernel_regularizer=regularizer_l2(1e-6),
              use_bias = T)
#################################################################################
# model <- keras_model_sequential() %>%
#   layer_dense(units = 128, activation = "relu", input_shape = c(dim(x.train)[2])
#               ,kernel_regularizer = regularizer_l2(0.01), 
#               bias_regularizer = regularizer_l2(0.01)
#               ) %>%
#   layer_dropout(rate = 0.2) %>%
#   layer_dense(units = 128, activation = "relu", use_bias = T              
#               ,kernel_regularizer = regularizer_l2(0.01), 
#               bias_regularizer = regularizer_l2(0.01)
#   ) %>%
#   layer_dropout(rate = 0.2) %>%
#   layer_dense(units = 64, activation = "relu", use_bias = T
#               ,kernel_regularizer = regularizer_l2(0.01), 
#               bias_regularizer = regularizer_l2(0.01)
#   ) %>%
#   layer_dropout(rate = 0.2) %>%
#   layer_dense(units = 32, activation = "relu", use_bias = T
#               ,kernel_regularizer = regularizer_l2(0.01), 
#               bias_regularizer = regularizer_l2(0.01)
#   ) %>%
#   layer_dropout(rate = 0.2) %>%
#   layer_dense(units = 16, activation = "relu", use_bias = T
#               ,kernel_regularizer = regularizer_l2(0.01), 
#               bias_regularizer = regularizer_l2(0.01)
#   ) %>%
#   layer_dropout(rate = 0.2) %>%
#   layer_dense(units = 10, activation = "relu", use_bias = T) %>%
#   layer_dense(units = dim(y.train)[2], activation = "tanh",
#               kernel_initializer=initializer_orthogonal(),
#               kernel_regularizer = regularizer_l2(0.01), use_bias = F) #model
#################################################################################
# Compile the model
model %>% compile(
  loss = "mse",
  optimizer = optimizer_rmsprop(learning_rate = 0.01),
  metrics = c("mse", "mae", "accuracy")
)
# summary(model)
# Train the model
history <- model %>% fit(
  x = x.train%>%as.matrix(), 
  y = y.train%>%as.matrix(),
  epochs = 50,
  batch_size = 64, 
  validation_split = 0.01
  # validation_data = list(x.val%>%as.matrix(), y.val%>%as.matrix())
)
# model %>% evaluate(x.test%>%as.matrix(), y.test%>%as.matrix())
# save_model_tf(model, "data/model-derivatives/model")
save_model_tf(model, "data/model-derivatives/model-lincs")
# model <- load_model_tf("data/model-derivatives/model")

# option 1
# mni <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/refs/mni_icbm152_nlin_sym_09c_CerebrA_nifti/annotated-xyz.rds")
# mni.ohe.all <-model.matrix(~ region - 1, 
#                        mni%>%mutate(region = as.factor(region)))
# colnames(mni.ohe.all) <- sub("region", "", colnames(mni.ohe.all))
# mni.hit <- mni %>%
#   select(starts_with("mni"))
# predictions <- predict(model, cbind(mni.hit,mni.ohe.all)%>%as.matrix())
# write_rds(predictions, "data/model-derivatives/predictions-cmap.rds")

# option 2 to ignore region label
mni <- readNIfTI("/Dedicated/jmichaelson-wdata/msmuhammad/refs/mni_icbm152_nlin_sym_09c_CerebrA_nifti/mni_icbm152_CerebrA_tal_nlin_sym_09c.nii", reorient = F)
origin <- c(97,133,79)
mni.hit <- which(mni>0, arr.ind = T) %>%
  as.data.frame() %>%
  mutate(mni_x = dim1-origin[1]) %>%
  mutate(mni_y = dim2-origin[2]) %>%
  mutate(mni_z = dim3-origin[3])
mni.hit.labeles <- mni.hit
mni <- readNIfTI("/Dedicated/jmichaelson-wdata/msmuhammad/refs/fsl-data/standard/MNI152_T1_1mm_brain.nii.gz", reorient = F)
origin <- c(91,127,73)
mni.hit <- which(mni>0, arr.ind = T) %>%
  as.data.frame() %>%
  mutate(mni_x = -1*(dim1-origin[1])) %>%
  mutate(mni_y = dim2-origin[2]) %>%
  mutate(mni_z = dim3-origin[3])
mni.hit.whole <- mni.hit
write_rds(mni.hit.whole, "data/mni-whole-brain-hit.rds")
# mni.hit <- full_join(mni.hit.labeles%>%select(starts_with("mni")),
#                      mni.hit.whole%>%select(starts_with("mni"))) %>%
#   distinct(mni_x, mni_y, mni_z)


predictions <- predict(model, mni.hit%>%select(starts_with("mni"))%>%as.matrix())
# write_rds(predictions, "data/model-derivatives/predictions-cmap_no-region-label.rds", compress="gz")
# write_rds(predictions, "data/model-derivatives/predictions-lincs_no-region-label_annot-pos.rds", compress="gz")
write_rds(predictions, "data/model-derivatives/predictions-lincs_no-region-label_whole-pos.rds", compress="gz")

# drug <- 1
# data <- data.frame(mni_x = mni.hit[,1], mni_y = mni.hit[,2], mni_z = mni.hit[,3], 
#                    correlation = predictions[,drug]) %>%
#   mutate(trans_corr = (((correlation-min(predictions[,drug]))/(max(predictions[,drug])-min(predictions[,drug])))*2)-1)
# data <- inner_join(data, mni %>% 
#                     select(starts_with("mni"), ends_with("region")))
# write_rds(data, "data/model-derivatives/data.rds")
#################################################################################
#################################################################################
#################################################################################
#################################################################################
# try a glm model for predictions
glm.predictions <- data.frame(mni.hit) %>% select(starts_with("mni"))
for (i in 5:ncol(all)) {
  # i=12
  drug <- colnames(all)[i]
  df <- all %>%
    rename(corr=i) %>%
    select(c(1:4),corr) %>%
    rename(mni_x=mni_x_allen) %>%
    rename(mni_y=mni_y_allen) %>%
    rename(mni_z=mni_z_allen) 
  # model <- glm(data = df, corr ~ mni_x + mni_y + mni_z + region)
  model <- glm(data = df, corr ~ mni_x + mni_y + mni_z)
  # tmp <- data.frame(orig = df$corr, pred = predict(model))
  # tmp %>% 
  #   ggplot(aes(x=orig, y=pred)) +
  #   geom_smooth()
  # # evaluate(p=as.vector(tmp$pred), a=as.vector(tmp$orig))
  # hist(residuals(model),breaks = 50)
  # glm.predictions$med <- predict(model, glm.predictions[,c("mni_x", "mni_y", "mni_z", "region")])
  glm.predictions$med <- predict(model, glm.predictions[,c("mni_x", "mni_y", "mni_z")])
  colnames(glm.predictions)[i-1] <- drug
  # print(i)
}
# write_rds(glm.predictions, "data/glm-predictions-cmap.rds")
# write_rds(glm.predictions, "data/glm-predictions-cmap_no-region-label.rds", compress="gz")
write_rds(glm.predictions, "data/glm-predictions-lincs_no-region-label_annot-pos.rds", compress="gz")
# write_rds(glm.predictions, "data/glm-predictions-lincs_no-region-label_whole-pos.rds", compress="gz")

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
use_condaenv("EDL")
library(keras)
k_clear_session()
set.seed(123)
#################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/brain-drug-map"
setwd(project.dir)
#################################################################################
load("data/exp/exp-corr-sex-age-ethn-donor-PMI_w-annot.rda")
cmap.of.int <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/LINCS/cmap.of.int.rds")
genes <- colnames(cmap.of.int)
#################################################################################
xyz <- annot.ls[,paste0("mni_", c("x", "y", "z"))] %>% as.matrix()
x.train <- xyz
y.train <- corr.brain.exp.ls %>% select(any_of(genes)) %>% as.matrix()
dim(y.train)
#################################################################################
# model
rm(model)
gc()
k_clear_session()
model <- keras_model_sequential(input_shape = c(dim(xyz)[2])) %>%
  layer_dense(units=512,activation="relu",
              use_bias=T,
  ) %>%
  layer_dropout(0.2) %>%
  layer_dense(units=512,activation="relu",
              use_bias=T,
  ) %>%
  layer_dense(units=256,activation="relu",
              use_bias=T,
  ) %>%
  layer_dense(units=64,activation="relu",
              use_bias=T,
  ) %>%
  layer_dropout(0.33) %>%
  layer_dense(units=ncol(y.train),activation="tanh",
              kernel_regularizer=regularizer_l1(1e-4),
              #	kernel_constraint=constraint_nonneg(),
              use_bias=T)
#################################################################################
# Compile the model
model %>% compile(
  loss = "mse",
  optimizer = optimizer_rmsprop(learning_rate = 0.01),
  metrics = c("mse", "mae"),
  loss_weights=list(c(1))
)
# summary(model)
# Train the model
history <- model %>% fit(
  x = x.train, 
  y = y.train,
  epochs = 100,
  batch_size = 256, 
  validation_split = 0.01,
  shuffle =T
  # validation_data = list(x.val%>%as.matrix(), y.val%>%as.matrix())
)

mni <- readNIfTI("/Dedicated/jmichaelson-wdata/msmuhammad/refs/mni_icbm152_nlin_sym_09c_CerebrA_nifti/mni_icbm152_CerebrA_tal_nlin_sym_09c.nii", reorient = F)
mni.pos <- which(mni>0, arr.ind = T)
origin <- c(97,133,79)
mni.pos <- mni.pos %>%
  as.data.frame() %>%
  mutate(x=dim1-origin[1]) %>%
  mutate(y=dim2-origin[2]) %>%
  mutate(z=dim3-origin[3]) %>%
  select(x,y,z) %>% as.matrix()
predictions <- predict(model, mni.pos)
#################################################################################
#################################################################################
