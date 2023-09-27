################################################################################
#       compute correlations between drug maps and BrainMap behavior maps      #
################################################################################
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
library(oro.nifti); library(ggh4x)
#################################################################################
#################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/brain-drug-map"
setwd(project.dir)
#################################################################################
# get mni brain coordinates
mni <- readNIfTI("/Dedicated/jmichaelson-wdata/msmuhammad/refs/fsl-data/standard/MNI152_T1_1mm_brain_mask.nii.gz", 
                 reorient = F)
origin <- c(91,127,73)
mni.hit <- which(mni>0, arr.ind = T) %>%
  as.data.frame() %>%
  mutate(mni_x = -1*(dim1-origin[1])) %>%
  mutate(mni_y = dim2-origin[2]) %>%
  mutate(mni_z = dim3-origin[3])
#################################################################################
# read behavior maps in MNI
beh.meta <- data.frame(file = list.files("/Dedicated/jmichaelson-wdata/msmuhammad/data/BrainMap/behavioral-domain/MNI/itk-registered/", full.names = F)) %>%
  mutate(domain = sub("\\..*", "", file)) %>%
  mutate(sub_domain = str_replace_all(sub("^[^.]+\\.(.*?)", "", 
                                          sub("_Z.nii.gz", "", file)),
                                      pattern = "[^A-Za-z]+", replacement = "_")) %>%
  mutate(full_label = paste0(domain, "_", sub_domain))
###
beh.maps.nii <- list()
beh.maps.df <- data.frame(mni.hit)
for (i in 1:nrow(beh.meta)) {
  nii <- readNIfTI(paste0("/Dedicated/jmichaelson-wdata/msmuhammad/data/BrainMap/behavioral-domain/MNI/itk-registered/", beh.meta$file[i]))
  beh.maps.nii[[i]] <- nii
  dim <- which(abs(nii)>0, arr.ind = T) %>%
    as.data.frame() %>%
    mutate(mni_x = -1*(dim1-origin[1])) %>%
    mutate(mni_y = dim2-origin[2]) %>%
    mutate(mni_z = dim3-origin[3])
  dim2 <- dim %>% as.data.frame() %>% mutate(int = nii[as.matrix(dim[,1:3])])
  colnames(dim2)[7] <- beh.meta$full_label[i]
  beh.maps.df <- inner_join(beh.maps.df, dim2)
  gc()
}
gc()
pdssave(beh.maps.df, file = "/Dedicated/jmichaelson-wdata/msmuhammad/data/BrainMap/behavioral-domain/MNI/itk-registered-all-beh-xyz-maps.rds")
# beh.maps.df <- pdsload("/wdata/msmuhammad/data/BrainMap/behavioral-domain/MNI/itk-registered-all-beh-xyz-maps.rds")
#################################################################################
# read drug activity maps
drugs.meta <- data.frame(file = list.files("data/maps/model-082223/all/", full.names = F)) %>%
  mutate(name = sub(".nii.gz", "", file))
drugs.maps.df <- pdsload("data/all-drug-map-predicted-exp-whole-brain-082223.rds.pxz")
#################################################################################
# compute correlation between both dataframes
t <- inner_join(beh.maps.df, drugs.maps.df)
gc()
drug.corr <- corr.func(t%>%select(beh.meta$full_label), 
                        t%>%select(drugs.meta$name),
                        method = "spearman")
pdssave(drug.corr, file = "data/drug-correlations-w-BrainMap-mni-maps.rds")
# drug.corr <- pdsload("data/drug-correlations-w-BrainMap-mni-maps.rds.pxz")
#################################################################################
# compute correlation between beh maps themselves 
beh.corr <- corr.func(t%>%select(beh.meta$full_label), 
                      t%>%select(beh.meta$full_label),
                      method = "spearman")
pdssave(beh.corr, file = "data/BrainMap-mni-maps-correlations-w-each-other.rds")
# beh.corr <- pdsload("data/BrainMap-mni-maps-correlations-w-each-other.rds.pxz")
#################################################################################
# compute correlation between drug maps
drug.act.corr <- corr.func(t%>%select(drugs.meta$name), 
                       t%>%select(drugs.meta$name),
                       method = "spearman")
pdssave(drug.act.corr, file = "data/drug-activity-maps-correlations-w-each-other.rds")
# drug.act.corr <- pdsload("data/drug-activity-maps-correlations-w-each-other.rds.pxz")
#################################################################################
inner_join(inner_join(beh.corr, beh.meta %>% select(V2 = full_label, 
                                                    domain_2 = domain, 
                                                    sub_domain_2 = sub_domain)), 
           beh.meta %>% rename(V1 = full_label)) %>% 
  filter(V1 != V2) %>%
  ggplot(aes(x=sub_domain, y=sub_domain_2, fill = r, label = ifelse(pval<0.01, "*", "")))+
  geom_tile()+
  geom_text(size = 3)+
  facet_grid2(cols = vars(reorder(domain, desc(domain))), rows = vars(domain_2), scales = "free", space = "free") +
  redblu.col.gradient + my.guides +
  labs(x="", y="", 
       title = "correlation between BrainMap behavior maps with each other",
       caption = "the maps were registered to MNI-space using itk-snap and filtered to keep positions mapped to the MNI brain mask") 



inner_join(drug.corr, beh.meta %>% rename(V1 = full_label)) %>%
  filter(V2 %in% c("methylphenidate", "sertraline", "venlafaxine", "fluoxetine",
                   "caffeine", "bupropion", "trazodone", "zolpidem",
                   "atomoxetine", "citalopram", "clonidine", "clonazepam",
                   "clozapine", "diazepam", "escitalopram", "fluvoxamine",
                   "haloperidol", "histamine", "ibuprofen", "ketamine",
                   "lidocaine", "nicotine", "orlistat", "paroxetine",
                   "risperidone", "sumatriptan", "topiramate")) %>%
  ggplot(aes(x=V2, y = sub_domain, fill = r, label = ifelse(pval<0.01, "*", "")))+
  geom_tile()+
  geom_text(size = 3)+
  facet_grid2(rows = vars(domain), scales = "free", space = "free") +
  redblu.col.gradient + my.guides +
  labs(x="", y="", 
       title = "correlation between predicted drug activity and BrainMap behavior maps",
       caption = paste0("the predicted drug activity per drug is based on the correlation between the transcriptomic signature per drug and the predicted gene expression profile from the DL model",
       "\n\ti.e., negative activity value means the drug is expected to make a change in that region",
       "\n\nthe BrainMap behavior maps were originally in Talairach space. They were then registered to the MNI-space using itk-snap.",
       "\n\nthe correlations here are based on the Spearman correlation between the predicted drug activity for a drug and the behavior MNI-registered map",
       "\n\ti.e., a negative correlation means we predicted the drug to make a change in a region that is activated for that specific behavior",
       "\n\ti.e., a positive correlation means we predicted the drug to have a similar gene expression to the region that is activated for that specific behavior"))





