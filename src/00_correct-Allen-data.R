################################################################################
#                  correct gene expression of Allen data                       #
################################################################################
rm(list = ls())
gc()
library(tidyverse)
set.seed(123)
#################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/brain-drug-map"
setwd(project.dir)
#################################################################################
#################################################################################
avg.brain.exp.ls <- read_rds("data/exp/combined-aggr-exp-of-all-subjects.rds") %>% as.data.frame()
annot.ls <- read_rds("data/exp/combined-annot-of-all-subjects.rds") %>% as.data.frame()
#################################################################################
# add the atlas labels to of the allen annotation
atlas.labels <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/refs/mni_icbm152_nlin_sym_09c_CerebrA_nifti/annotated-xyz.rds")
annot.ls <- left_join(annot.ls %>%
                   rename(mni_x_allen=mni_x)  %>%
                   rename(mni_y_allen=mni_y)  %>%
                   rename(mni_z_allen=mni_z)  %>%
                   mutate(mni_x = round(mni_x_allen))  %>%
                   mutate(mni_y = round(mni_y_allen))  %>%
                   mutate(mni_z = round(mni_z_allen)),
                 atlas.labels)
# tmp.2 <- cbind(tmp, avg.brain.exp.ls) %>% drop_na(region)
# avg.brain.exp.ls.filt <- tmp.2[,29:ncol(tmp.2)]
# annot.ls.filt <- tmp.2[,1:28]
#################################################################################
# check if there's a correlation between measured gene exp and sex, age, ethnicity, PMI
summary(lm(data = data.frame(exp = avg.brain.exp.ls[,100], annot.ls %>% select(donor, sex, age, ethnicity, PMI)),
           exp ~ sex)) #sig for some genes
summary(lm(data = data.frame(exp = avg.brain.exp.ls[,1], annot.ls %>% select(donor, sex, age, ethnicity, PMI)),
           exp ~ age)) #sig
summary(lm(data = data.frame(exp = avg.brain.exp.ls[,1], annot.ls %>% select(donor, sex, age, ethnicity, PMI)),
           exp ~ PMI)) #sig
summary(lm(data = data.frame(exp = avg.brain.exp.ls[,1], annot.ls %>% select(donor, sex, age, ethnicity, PMI)),
           exp ~ ethnicity)) #sig
summary(lm(data = data.frame(exp = avg.brain.exp.ls[,1], annot.ls %>% select(donor, sex, age, ethnicity, PMI)),
           exp ~ donor)) #sig
summary(lm(data = data.frame(exp = avg.brain.exp.ls[,1], annot.ls %>% select(donor, sex, age, ethnicity, PMI)),
           exp ~ sex+age+PMI+donor+ethnicity))
# from these lm results, you should correct for these
#################################################################################
# correct avg.exp for each gene based on these covariates
library(doMC)
registerDoMC(cores = 8)
corr.brain.exp.ls <- foreach(i = 1:ncol(avg.brain.exp.ls), .combine = cbind) %dopar% {
  # i=12
  gene <- colnames(avg.brain.exp.ls)[i]
  df <- data.frame(exp = avg.brain.exp.ls[,i], 
                   annot.ls) %>%
    mutate(sc_exp = scale(exp, scale = T, center = T)[,1])
  df <- df%>%
    mutate(corr_exp = residuals(lm(data = df, sc_exp ~ sex+age+PMI+donor+ethnicity)))
  df.2 <- left_join(annot.ls, df) %>%
    select(corr_exp)
  colnames(df.2) <- gene
  return(df.2)
}
save(corr.brain.exp.ls, avg.brain.exp.ls, annot.ls, file = "data/exp/exp-corr-sex-age-ethn-donor-PMI_w-annot.rda")
#################################################################################