################################################################################
#        compute correlations between drug maps and neurosynth features        #
################################################################################
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
#################################################################################
#################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/brain-drug-map"
setwd(project.dir)
#################################################################################
# extras
brain.networks <- c("default mode", "sensorimotor", "visual", "limbic", "dmn", 
                    "central executive", "salience", "dorsal attention",
                    "fronto parietal", "language", "cingulo opercular"
                    ,"network"
)
brain.regions <- c("putamen", "cortex", "basal ganglia", "caudate", "thalamus", "insula", "cerebellar",
                   "striatum", "promotor cortex", "frontal cortices", "prefrontal cortex", "amygdala",
                   "cingulate", "anterior cingulate", "cingulate gyrus")
drugs <- c("methylphenidate", "sertraline", "venlafaxine", "fluoxetine",
           "caffeine", "bupropion", "trazodone", "zolpidem",
           "atomoxetine", "citalopram", "clonidine", "clonazepam",
           "clozapine", "diazepam", "escitalopram", "fluvoxamine",
           "haloperidol", "histamine", "ibuprofen", "ketamine",
           "lidocaine", "nicotine", "orlistat", "paroxetine",
           "risperidone", "sumatriptan", "topiramate")
concepts <- read_csv("/Dedicated/jmichaelson-wdata/msmuhammad/data/cognitive-atlas/concepts.csv")
phenotypes <- read_csv("/Dedicated/jmichaelson-wdata/msmuhammad/data/cognitive-atlas/phenotypes.csv")
#################################################################################
# load neurosynth MNI features
ns.features.all.0 <- pdsload(fname = "data/neurosynth/features-data-by-xyz-in-MNI-only.rds.pxz") %>% as.data.frame()
# load drug maps
drug.corr <- pdsload(fname = "data/all-drug-map-predicted-exp-whole-brain-082223.rds.pxz") %>% as.data.frame()
gc()
#################################################################################
# keep features of interest
# keep the manualy selected ones
tmp <- read_csv("data/neurosynth/features-dictionary.csv") %>% 
  filter(beh==1)
ns.features.all <- ns.features.all.0 %>%
  select(x, y, z, any_of(c(brain.networks, tmp$f)))
ns.features <- ns.features.all # keep them all
# ns.features <- ns.features.all[,c(T,T,T,grepl(str_c(tolower(brain.networks), collapse = "|"), 
#                                               colnames(ns.features.all))[-c(1:3)])]

# get voxels in common between drug correlations and neurosynth feature space
dimensions <- right_join(ns.features %>% select(1:3),
                         drug.corr %>% select(x = mni_x, y=mni_y, z=mni_z)) %>%
  distinct()
m1 <- inner_join(dimensions, drug.corr %>% 
                   select(-c(starts_with("dim"))) %>% 
                   rename(x = mni_x, y = mni_y, z = mni_z)) %>% 
  select(-c(x,y,z))
# missing positions
miss <- left_join(anti_join(dimensions, ns.features %>% 
                              mutate_at(.vars = vars(x,y,z), .funs = function(x) round(x)) %>%
                              distinct(x,y,z, .keep_all = T)),
                  ns.features) %>% 
  mutate_all(.funs = function(x) replace_na(x, 0))
gc()
# select features of interest
m2 <- rbind(inner_join(dimensions, ns.features %>% 
                         mutate_at(.vars = vars(x,y,z), .funs = function(x) round(x)) %>%
                         distinct(x,y,z, .keep_all = T)),
            miss) %>%
  select(-c(x,y,z))
gc()
drug.f.maps <- corr.func(m1, m2, method = "spearman")
gc()
# pdssave(drug.f.maps, file = "data/drug-features-correlations-by-xyz-in-predicted-MNI-wpval.rds")
pdssave(drug.f.maps, file = "data/drug-features-correlations-by-xyz-in-predicted-MNI-wpval-beh.rds")
# drug.f.maps <- pdsload(fname = "data/drug-features-correlations-by-xyz-in-predicted-MNI.rds.pxz")
#################################################################################
gc()
long.drug.f.map <- drug.f.maps %>%
  as.data.frame() %>%
  rownames_to_column("c_drug") %>%
  pivot_longer(cols = colnames(ns.features)[4:ncol(ns.features)], names_to = "feature", values_to = "val")

long.drug.f.map %>%
  filter(c_drug %in% drugs ,
         # grepl(brain.networks, feature)
         # grepl(tolower(brain.regions), tolower(feature))
         # grepl(str_c(tolower(brain.networks), collapse = "|"), feature)
         grepl(str_c(tolower(brain.regions), collapse = "|"), feature)
         # grepl(str_c(tolower(phenotypes$phenotype[which(phenotypes$category == "disorder")]), collapse = "|"), feature)
         # concepts categories are: action, attention, emotion, executive-cognitive control, language
         # memory, motivation, perception, reasoning and decision making, social function
         # grepl(str_c(tolower(concepts$concept[which(concepts$category == "motivation")]), collapse = "|"), feature)
         ) %>%
  ggplot(aes(x=c_drug, y=feature, fill = val)) +
  geom_tile()+
  # geom_text(size = 2) +
  redblu.col.gradient +
  null_labs +
  my.guides
write_csv(long.drug.f.map %>% filter(c_drug == "methylphenidate"), file = "data/MPH-features-correlations-V2.csv")

#################################################################################

