################################################################################
#              predeicting activity of drugs in brain xyz positions            #
################################################################################
rm(list = ls())
gc()
# source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
library(tidyverse)
set.seed(123)
#################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/brain-drug-map"
setwd(project.dir)
#################################################################################
load("data/exp/exp-corr-sex-age-ethn-PMI_w-annot.rda")
gc()
atlas.annot <- read_rds("../../refs/mni_icbm152_nlin_sym_09c_CerebrA_nifti/annotated-xyz.rds")
atlas.mni <- read_rds("../../refs/mni_icbm152_nl_VI_nifti/brain-pos.rds")
#################################################################################
# drug.sig <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/LINCS/LINCS-combined-corrected-14.rds") %>%
#   dplyr::filter(symbol %in% colnames(corr.brain.exp.ls)) %>%
#   column_to_rownames("symbol")
drug.sig <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/LINCS/cmap.of.int_v2.rds") %>%
  t()
table(rownames(drug.sig) %in% colnames(corr.brain.exp.ls))
genes <- rownames(drug.sig)[(rownames(drug.sig)%in%colnames(corr.brain.exp.ls))]

scaled.brain.exp.2 <- scale(corr.brain.exp.ls,scale = T, center = T)
scaled.brain.exp.2 <- scaled.brain.exp.2[,genes] %>% t()
drug.sig <- drug.sig[genes,]
all(rownames(scaled.brain.exp.2) == rownames(drug.sig))
# ignore scaling gene expression in drugs, it actually makes sense and a difference
drug.sig <- t(scale(t(drug.sig), scale = T, center = T))
dim(scaled.brain.exp.2)
dim(drug.sig)
#################################################################################
# compute correlations with drugs
# I transform the drug correlations
drug.corr <- cor(scaled.brain.exp.2, drug.sig)
drug.corr.annot <- cbind(annot.ls, drug.corr)

# you don't need scaling?
# scaled.drug.corr <- scale(drug.corr, scale = T, center = T)
trans.drug.corr <- drug.corr
for (i in 1:ncol(drug.corr)) {
  d.cor <- drug.corr[,i]
  min.scaled <- min(d.cor)
  max.scaled <- max(d.cor)
  trans.drug.corr[,i] <- (d.cor - min.scaled) / (max.scaled - min.scaled) * 2 - 1
}
trans.drug.corr.annot <- cbind(annot.ls, trans.drug.corr)
#################################################################################
# you now have the data needed for any model
# you want 3 models with different approaches
# Model 1: takes drug correlations, allen xyz positions, corresponding region annot from atlas
#           predicts drug corr for all regions with xyz from region annot atlas
# Model 2: takes drug correlations, allen xyz positions
#           predicts drug corr for all xyz positions from whole MNI brain
# Model 3: takes drug correlations, allen xyz positions, corresponding region annot from atlas and add "NOT" region for unannotated regions from allen
#           predicts drug corr for all regions with xyz from whole MNI brain with named regions from atlas annot
#################################################################################
################################ Model 1 ########################################
#################################################################################
m1.trans.drug.corr.annot <- trans.drug.corr.annot %>%
  drop_na(region)
m1.xyz <- m1.trans.drug.corr.annot %>%
  select(mni_x, mni_y, mni_z, region)
m1.drug.corr <- m1.trans.drug.corr.annot %>%
  select(colnames(drug.sig))
# try a glm model for predictions
m1.glm.predictions <- data.frame(atlas.annot) %>% select(starts_with("mni"), region, h_region)
for (i in (ncol(annot.ls)+1):ncol(m1.trans.drug.corr.annot)) {
  # i=29
  drug <- colnames(m1.trans.drug.corr.annot)[i]
  df <- m1.trans.drug.corr.annot %>%
    rename(corr=i) %>%
    select(mni_x=mni_x_allen, mni_y=mni_y_allen, mni_z=mni_z_allen,region,corr) %>%
    mutate(region=factor(region))
  model <- glm(data = df, corr ~ mni_x + mni_y + mni_z + region)
  
  m1.glm.predictions$med <- predict(model, m1.glm.predictions[,c("mni_x", "mni_y", "mni_z", "region")])
  colnames(m1.glm.predictions)[ncol(m1.glm.predictions)] <- drug
  # print(i)
}
write_rds(m1.glm.predictions, "data/model-derivatives/01_m1-glm-lincs-predictions.rds")
#################################################################################
################################ Model 2 ########################################
#################################################################################
m2.trans.drug.corr.annot <- trans.drug.corr.annot
m2.xyz <- m2.trans.drug.corr.annot %>%
  select(mni_x, mni_y, mni_z)
m2.drug.corr <- m2.trans.drug.corr.annot %>%
  select(colnames(drug.sig))
# try a glm model for predictions
m2.glm.predictions <- left_join(data.frame(atlas.mni) %>% select(starts_with("mni")),
                                atlas.annot %>% select(starts_with("mni"), ends_with("region")))
for (i in (ncol(annot.ls)+1):ncol(m2.trans.drug.corr.annot)) {
  # i=29
  drug <- colnames(m2.trans.drug.corr.annot)[i]
  df <- m2.trans.drug.corr.annot %>%
    rename(corr=i) %>%
    select(mni_x=mni_x_allen, mni_y=mni_y_allen, mni_z=mni_z_allen,corr)
  model <- glm(data = df, corr ~ mni_x + mni_y + mni_z)
  
  m2.glm.predictions$med <- predict(model, m2.glm.predictions[,c("mni_x", "mni_y", "mni_z")])
  colnames(m2.glm.predictions)[ncol(m2.glm.predictions)] <- drug
  # print(i)
}
write_rds(m2.glm.predictions, "data/model-derivatives/01_m2-glm-lincs-predictions.rds")
#################################################################################
################################ Model 3 ########################################
#################################################################################
m3.trans.drug.corr.annot <- trans.drug.corr.annot %>%
  mutate(region=ifelse(is.na(region), "NOT", region)) %>%
  mutate(h_region=ifelse(is.na(h_region), "NOT", h_region))
m3.xyz <- m3.trans.drug.corr.annot %>%
  select(mni_x, mni_y, mni_z, region)
m3.drug.corr <- m3.trans.drug.corr.annot %>%
  select(colnames(drug.sig))
# try a glm model for predictions
m3.glm.predictions <- full_join(data.frame(atlas.mni) %>% select(starts_with("mni")),
                                atlas.annot %>% select(starts_with("mni"), ends_with("region"))) %>%
  mutate(region=ifelse(is.na(region), "NOT", region)) %>%
  mutate(h_region=ifelse(is.na(h_region), "NOT", h_region))
for (i in (ncol(annot.ls)+1):ncol(m3.trans.drug.corr.annot)) {
  # i=29
  drug <- colnames(m3.trans.drug.corr.annot)[i]
  df <- m3.trans.drug.corr.annot %>%
    rename(corr=i) %>%
    select(mni_x=mni_x_allen, mni_y=mni_y_allen, mni_z=mni_z_allen,region,corr) %>%
    mutate(region=factor(region))
  model <- glm(data = df, corr ~ mni_x + mni_y + mni_z + region)
  
  m3.glm.predictions$med <- predict(model, m3.glm.predictions[,c("mni_x", "mni_y", "mni_z", "region")])
  colnames(m3.glm.predictions)[ncol(m3.glm.predictions)] <- drug
  # print(i)
}
write_rds(m3.glm.predictions, "data/model-derivatives/01_m3-glm-lincs-predictions.rds")
#################################################################################
m1.m2.meph <- inner_join(m1.glm.predictions %>% 
                           select(starts_with("mni"), ends_with("region"),m1_meph=methylphenidate),
                         m2.glm.predictions %>% 
                           select(starts_with("mni"), ends_with("region"),m2_meph=methylphenidate))
m1.m3.meph <- inner_join(m1.glm.predictions %>% 
                           select(starts_with("mni"), ends_with("region"),m1_meph=methylphenidate),
                         m3.glm.predictions %>% 
                           select(starts_with("mni"), ends_with("region"),m3_meph=methylphenidate))
m2.m3.meph <- inner_join(m2.glm.predictions %>% 
                           select(starts_with("mni"), ends_with("region"),m2_meph=methylphenidate),
                         m3.glm.predictions %>% 
                           select(starts_with("mni"), ends_with("region"),m3_meph=methylphenidate))
p12 <- m1.m2.meph %>%
  ggplot(aes(x=m1_meph, y=m2_meph)) +
  geom_smooth()
p13 <- m1.m3.meph %>%
  ggplot(aes(x=m1_meph, y=m3_meph)) +
  geom_smooth()
p23 <- m2.m3.meph %>%
  ggplot(aes(x=m2_meph, y=m3_meph)) +
  geom_smooth()
cowplot::plot_grid(p12, p13, p23)
#################################################################################
#################################################################################
############################## Deep Learning ####################################
#################################################################################
#################################################################################
# model 4 is the same as model 2, but with a DL model
library(tensorflow)
use_condaenv("EDL")
library(keras)
set.seed(123)
# Define the neural network architecture
rm(model)
k_clear_session()

m4.trans.drug.corr.annot <- trans.drug.corr.annot
m4.xtrain <- m4.trans.drug.corr.annot %>%
  select(mni_x=mni_x_allen, mni_y=mni_y_allen, mni_z=mni_z_allen)
m4.ytrain <- m4.trans.drug.corr.annot %>%
  select(colnames(drug.sig))

model <- keras_model_sequential(input_shape = c(dim(m4.xtrain)[2])) %>%
  layer_dense(units = 256, activation = "relu") %>%
  layer_dense(units = 256, activation = "relu") %>%
  layer_dense(units = 64, activation = "relu",
              use_bias = TRUE,
              kernel_regularizer = regularizer_l2(1e-6),
  ) %>%
  layer_dense(units = ncol(m4.ytrain), activation = "tanh", 
              kernel_regularizer=regularizer_l2(1e-6),
              use_bias = T)
# Compile the model
model %>% compile(
  loss = "mse",
  optimizer = optimizer_rmsprop(learning_rate = 0.01),
  metrics = c("mse", "mae")
)
# summary(model)
# Train the model
history <- model %>% fit(
  x = m4.xtrain%>%as.matrix(), 
  y = m4.ytrain%>%as.matrix(),
  epochs = 20,
  batch_size = 64, 
  validation_split = 0.01
  # validation_data = list(x.val%>%as.matrix(), y.val%>%as.matrix())
)
save_model_tf(model, "data/model-derivatives/model-cmap")

m4.xpred <- atlas.mni %>%
  select(starts_with("mni")) %>%
  as.matrix()
m4.dl.predictions <- predict(model, m4.xpred)
colnames(m4.dl.predictions) <- colnames(drug.sig)
m4.dl.predictions <- cbind(left_join(atlas.mni, atlas.annot%>%rename(annot_dim1=dim1,annot_dim2=dim2,annot_dim3=dim3)), m4.dl.predictions)
write_rds(m4.dl.predictions, "data/model-derivatives/01_m4-dl-cmap-predictions.rds")
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
