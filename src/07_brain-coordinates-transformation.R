################################################################################
#               convert brain MNI positions to tal using GngerALE              #
################################################################################
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
library(oro.nifti)
#################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/brain-drug-map"
setwd(project.dir)
#################################################################################
# tried converting tal to mni
# had an issue in finding a way to only extract brain?
tal <- readNIfTI("/wdata/msmuhammad/data/BrainMap/behavioral-domain/Tal/Action.Execution.All_Z.nii.gz", 
                 reorient = F)
dim(tal)
origin.tal <- c(41,59,32)
# not correct. you have intensities with negative values!
tal.hit <- which(abs(tal)>0, arr.ind = T) %>%
  as.data.frame() %>%
  mutate(tal_x = -1*(dim1-origin.tal[1])) %>%
  mutate(tal_y = dim2-origin.tal[2]) %>%
  mutate(tal_z = dim3-origin.tal[3])

####

# try getting mni xyz and save it for conversion by GingerALE
mni <- readNIfTI("/Dedicated/jmichaelson-wdata/msmuhammad/refs/fsl-data/standard/MNI152_T1_1mm_brain_mask.nii.gz", 
                 reorient = F)
origin <- c(91,127,73)
mni.hit <- which(mni>0, arr.ind = T) %>%
  as.data.frame() %>%
  mutate(mni_x = -1*(dim1-origin[1])) %>%
  mutate(mni_y = dim2-origin[2]) %>%
  mutate(mni_z = dim3-origin[3])
# save mni positions in a txt file to be used by GingerALE GUI
write.table(mni.hit[,4:6], row.names = F, col.names = F, "data/mni-brain-coordinates.txt")

####

# read the converted coordinates and round them up
mni2tal <- read_table("data/mni-brain-coordinates-in-tal-space-by-GingerALE.txt", 
                      col_names = c("Tal_x_r", "Tal_y_r", "Tal_z_r"))
mni.tal.full <- cbind(mni.hit, mni2tal,t) %>%
  mutate(tal_x = round(Tal_x_r),
         tal_y = round(Tal_y_r),
         tal_z = round(Tal_z_r))

####

# read mango tranformation mat
mango.mni2tal.mat <- read_table("/wdata/msmuhammad/data/BrainMap/behavioral-domain/transformation-mat/mni2tal-mango.mat", col_names = F)
mango.tal2mni.mat <- read_table("/wdata/msmuhammad/data/BrainMap/behavioral-domain/transformation-mat/tal2mni-mango.mat", col_names = F)

talmni.mango <- as.matrix(tal.hit[,1:4])%*%as.matrix(mango.tal2mni.mat)
mnital.mango <- as.matrix(cbind(mni.hit[,4:6],1))%*%as.matrix(mango.mni2tal.mat)
