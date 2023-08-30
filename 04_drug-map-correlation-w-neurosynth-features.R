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
#################################################################################
# load neurosynth MNI features
ns.features <- pdsload(fname = "data/neurosynth/features-data-by-xyz-in-MNI-only.rds.pxz") %>% as.data.frame()
# load drug maps
drug.corr <- pdsload(fname = "data/all-drug-map-predicted-exp-whole-brain-082223.rds.pxz") %>% as.data.frame()
#################################################################################
# get voxels in common between drug correlations and neurosynth feature space
dimensions <- inner_join(ns.features %>% select(1:3),
                         drug.corr %>% select(x = mni_x, y=mni_y, z=mni_z)) %>%
  distinct()
drug.f.maps <- cor(inner_join(dimensions, drug.corr %>% 
                                select(-c(starts_with("dim"))) %>% 
                                rename(x = mni_x, y = mni_y, z = mni_z)) %>% 
                     select(-c(x,y,z)),
                   inner_join(dimensions, ns.features %>% 
                                mutate_at(.vars = vars(x,y,z), .funs = function(x) round(x)) %>%
                                distinct(x,y,z, .keep_all = T)) %>%
                     select(-c(x,y,z)),
                   method = "spearman")
pdssave(drug.f.maps, file = "data/drug-features-correlations-by-xyz-in-MNI-only.rds")
#################################################################################
long.drug.f.map <- drug.f.maps %>%
  as.data.frame() %>%
  rownames_to_column("c_drug") %>%
  pivot_longer(cols = colnames(ns.features)[4:ncol(ns.features)], names_to = "feature", values_to = "val")

long.drug.f.map %>%
  filter(c_drug %in% drugs, 
         grepl(str_c(tolower(brain.networks), collapse = "|"), feature)
         # grepl(str_c(tolower(personality.traits), collapse = "|"), feature)
         # grepl(str_c(tolower(brain.regions), collapse = "|"), feature)
         ) %>%
  ggplot(aes(x=c_drug, y=feature, fill = val, label = round(val, 3))) +
  geom_tile()+
  geom_text(size = 2) +
  redblu.col.gradient +
  null_labs +
  my.guides
write_csv(long.drug.f.map %>% filter(c_drug == "methylphenidate"), file = "data/MPH-features-correlations.csv")

#################################################################################
# extras
brain.networks <- c("default mode", "sensorimotor", "visual", "limbic", 
                    "central executive", "salience", "dorsal attention",
                    "fronto parietal", "language", "cingulo opercular"
                    # ,"network"
                    )
personality.traits <- c("agreeableness","alexithymia","conscientiousness",
                        "cooperativeness","extraversion","harm avoidance",
                        "honesty","impulsiveness","neuroticism",
                        "novelty seeking","openness to experience",
                        "persistence","reward dependence",
                        "self directedness","self esteem","self transcendence",
                        "sensation-seeking")
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
