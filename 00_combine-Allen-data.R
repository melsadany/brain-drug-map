################################################################################
#      building glm per gene per tissue for all snps and extract weights       #
################################################################################
rm(list = ls())
gc()
# source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
library(tidyverse)
#################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/brain-drug-map"
setwd(project.dir)

#################################################################################
donors <- c(9861, 10021, 12876, 14380, 15496, 15697)
# mapp the nii file names to donor names in microarray
# they take the same order of the line above
nii.donors <- c("H0351.2001", "H0351.2002", "H0351.1009", "H0351.1012", "H0351.1015", "H0351.1016")
donors.meta <- data.frame(donor = donors, nii_donor = nii.donors, sex =c("M", "M", "M", "M", "F", "M"),
                          age = c(24, 39, 57, 31, 49, 55), #years
                          ethnicity = c("BAA", "BAA", "WC", "WC", "H", "WC"), # WC white/Caucasian BAA Black/African-American H Hispanic
                          PMI = c(23, 10, 26, 17, 30, 18)) #hours
write_tsv(donors.meta, "data/donors-metadata.tsv")

# only done once
avg.brain.exp.ls <- list()
annot.ls <- list()
for (i in 1:length(donors)) {
  # i=1
  donor <- donors[i]
  nii.donor <- nii.donors[i]
  # annot well-id are unique across samples and donors
  annot <- read_csv(paste0("data/exp/normalized_microarray_donor", donor, "/SampleAnnot.csv"))
  # probes are the same across all samples and donors
  probes <- read_csv(paste0("data/exp/normalized_microarray_donor", donor, "/Probes.csv"))
  brain.exp <- read_csv(paste0("data/exp/normalized_microarray_donor", donor, "/MicroarrayExpression.csv"), col_names = c("probe", annot$well_id))
  tmp <- cbind(gene = probes$gene_symbol[which(probes$probe_id %in% brain.exp$probe)], brain.exp[,-1])
  avg.brain.exp <- aggregate(tmp[,-1], list(tmp$gene), FUN = mean)
  avg.brain.exp.2 <- avg.brain.exp %>%
    column_to_rownames("Group.1") %>%
    t()
  avg.brain.exp.ls[[as.character(donor)]] <- avg.brain.exp.2
  
  # # read the nifti image, get the origin, transform the origin to 000 xyz position, and transform all the mri positions by adding the difference
  # library(ANTsRCore, lib.loc = "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/miniconda3/envs/ENA/lib/R/library")
  # library(ANTsR, lib.loc = "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/miniconda3/envs/ENA/lib/R/library")
  # donor.t1 <- antsImageRead(paste0("/Dedicated/jmichaelson-wdata/msmuhammad/projects/Allen/human_brain/", nii.donor, "/sub-", nii.donor, "_T1.nii.gz"))
  # ants.origin <- antsGetOrigin(donor.t1)
  # t1.origins <- read_csv("/Dedicated/jmichaelson-wdata/msmuhammad/data/Allen/Imaging/T1-origins-fsleyes.csv")
  # origin <- t1.origins %>%
  #   filter(subject == paste0("Sub-", nii.donor)) %>%
  #   select(x,y,z) %>%
  #   unlist() 
  # mni.origin <- t1.origins %>%
  #   filter(subject == "Mni-allen") %>%
  #   select(x,y,z) %>%
  #   unlist() 
  # 
  # hit <- which(donor.t1>0,arr.ind=T)
  # # I'm following this equation
  # # MRI+MNI = diff
  # diff <- c(origin[[1]]-mni.origin[[1]], origin[[2]]-mni.origin[[2]], origin[[3]]-mni.origin[[3]])
  # 
  # annot.new <- annot %>%
  #   mutate(mni_x_calc = -(diff[1]-mri_voxel_x)) %>%
  #   mutate(mni_y_calc = -(diff[2]-mri_voxel_y)) %>%
  #   mutate(mni_z_calc = -(diff[3]-mri_voxel_z)) %>%
  #   mutate(mni_x_ants_calc = (ants.origin[1]+mri_voxel_x)) %>%
  #   mutate(mni_y_ants_calc = (ants.origin[2]+mri_voxel_y)) %>%
  #   mutate(mni_z_ants_calc = (ants.origin[3]+mri_voxel_z)) %>%
  #   mutate(mni_x_0_calc = mri_voxel_x-origin[[1]]) %>%
  #   mutate(mni_y_0_calc = mri_voxel_y-origin[[2]]) %>%
  #   mutate(mni_z_0_calc = mri_voxel_z-origin[[3]])
  
  annot.ls[[as.character(donor)]] <- left_join(annot %>%
                                                 mutate(donor = donor), 
                                               donors.meta %>% 
                                                 filter(donor==donor))
}

avg.brain.exp.ls <- do.call(rbind, avg.brain.exp.ls)
annot.ls <- do.call(rbind, annot.ls)

write_rds(avg.brain.exp.ls, "data/exp/combined-aggr-exp-of-all-subjects.rds")
# avg.brain.exp.ls <- read_rds("data/exp/combined-aggr-exp-of-all-subjects.rds") %>% as.data.frame()
write_rds(annot.ls, "data/exp/combined-annot-of-all-subjects.rds")
# annot.ls <- read_rds("data/exp/combined-annot-of-all-subjects.rds") %>% as.data.frame()


#################################################################################
