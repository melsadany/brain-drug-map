################################################################################
#           produce interactive 3d maps for drug maps and save them            #
################################################################################
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
library(plotly)
library(ggseg3d)
library(ggsegDKT)
library(ggseg)
scene <- list(camera = list(eye = list(x = 1.5, y = 0, z = 0)))
#################################################################################
#################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/brain-drug-map"
setwd(project.dir)
#################################################################################
#################################################################################
# read the drug activity maps
annot <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/refs/mni_icbm152_nlin_sym_09c_CerebrA_nifti/annotated-xyz.rds")[,5:9]
drug.maps <- pdsload("data/all-drug-map-predicted-exp-whole-brain-082223.rds.pxz")
drugs <- colnames(drug.maps)[7:ncol(drug.maps)]
# #### 
# # indep section to save average correlations per region in a csv
# avg.drug.map <- inner_join(annot,
#                        drug.maps[,-c(1:3)]) %>%
#   as.data.frame() %>%
#   mutate(region_2 = sub("_", " ", region)) %>%
#   mutate(h_region = sub("rh_", "right ", h_region)) %>%
#   mutate(h_region = sub("lh_", "left ", h_region)) %>%
#   mutate(h_region = sub("_", " ", h_region)) %>%
#   pivot_longer(cols = drugs, names_to = "drug", values_to = "corr") %>%
#   group_by(drug, h_region) %>%
#   dplyr::summarize(mean_corr = mean(corr)) %>%
#   drop_na() 
# gc()
# tmp <- avg.drug.map %>% pivot_wider(names_from="drug", values_from="mean_corr")
# write_csv(tmp, file = "data/all-drug-map-activity-averaged-by-anatomical-region-082223.csv")
# ####
gc()
#################################################################################
# iterate over the list of drugs and make a map for each
registerDoMC(cores = 6)
foreach(i=1:length(drugs)) %dopar% {
  gc()
  med <- drugs[i]
  drug.map <- inner_join(annot,
                         cbind(drug.maps[,4:6], drug.maps[,med])) %>%
    as.data.frame() %>%
    rename(corr=6) %>%
    mutate(region_2 = sub("_", " ", region)) %>%
    mutate(h_region = sub("rh_", "right ", h_region)) %>%
    mutate(h_region = sub("lh_", "left ", h_region)) %>%
    mutate(h_region = sub("_", " ", h_region))
  gc()
  #################################################################################
  # averaging data based on h_region label
  avg.h.data <- drug.map %>%
    group_by(h_region) %>%
    dplyr::summarize(mean_corr = mean(corr)) %>%
    drop_na()
  gc()
  #################################################################################
  # left hemisphere
  data.3d.lh <- dkt_3d %>%
    filter(surf == "inflated" & hemi == "left") %>% 
    unnest(ggseg_3d) %>% 
    ungroup() %>% 
    select(region) %>% 
    na.omit()
  data.3d.lh <- left_join(data.3d.lh, 
                          avg.h.data %>%
                            filter(!grepl("right", h_region)) %>%
                            mutate(region = sub("left ", "", h_region)))
  gc()
  #################################################################################
  # right hemisphere
  data.3d.rh <- dkt_3d %>%
    filter(surf == "inflated" & hemi == "right") %>% 
    unnest(ggseg_3d) %>% 
    ungroup() %>% 
    select(region) %>% 
    na.omit()
  data.3d.rh <- left_join(data.3d.rh, 
                          avg.h.data %>%
                            filter(!grepl("left", h_region)) %>%
                            mutate(region = sub("right ", "", h_region)))
  gc()
  #################################################################################
  #  subcortical regions
  h.data.3d <- aseg_3d %>% 
    unnest(cols = ggseg_3d) %>% 
    select(label) %>%
    mutate(atlas_label=label) %>%
    mutate(label = tolower(label)) %>%
    mutate(label = sub("-", " ", label)) %>%
    mutate(label = sub("-", " ", label)) %>%
    mutate(label = sub("-", " ", label)) %>%
    mutate(label = sub("_", " ", label)) %>%
    mutate(label = sub("_", " ", label)) %>%
    arrange(label)
  
  # this is me editing the labels manually to match same annotation of the data that I have
  h.data.3d$region <- c("left third ventricle", "left fourth ventricle",
                        "left brainstem", rep(NA,5), h.data.3d$label[9:11],
                        "left cerebellum gray matter", h.data.3d$label[13:14],
                        "left inferior lateral ventricle", 
                        h.data.3d$label[16:18], "left thalamus",
                        "left ventral diencephalon",
                        h.data.3d$label[21:23],
                        "right cerebellum gray matter", h.data.3d$label[25:26],
                        "right inferior lateral ventricle", 
                        h.data.3d$label[28:30], "right thalamus",
                        "right ventral diencephalon")
  h.data.3d.2 <- left_join(h.data.3d, avg.h.data %>% rename(region = h_region)) %>%
    # filter(!grepl("Ventricle|Putamen|Amygdala", label)) %>% 
    drop_na() %>%
    select(-region, -label) %>%
    rename(label=atlas_label)
  gc()
  #################################################################################
  # combine plots
  lh <- ggseg3d(data.3d.lh, 
                atlas = dkt_3d,
                colour = "mean_corr", text = "mean_corr",
                palette = c("forestgreen" = min(avg.h.data$mean_corr), 
                            "white" = 0, 
                            "firebrick" = max(avg.h.data$mean_corr)), 
                hemisphere = "left") %>%
    remove_axes()
  rh <- ggseg3d(data.3d.rh, 
                atlas = dkt_3d,
                colour = "mean_corr", text = "mean_corr",
                palette = c("forestgreen" = min(avg.h.data$mean_corr), 
                            "white" = 0, 
                            "firebrick" = max(avg.h.data$mean_corr)), 
                hemisphere = "right",
                show.legend = F) %>%
    remove_axes()
  sc <- ggseg3d(h.data.3d.2, atlas = aseg_3d, 
                colour = "mean_corr", text = "mean_corr", 
                na.alpha= .5,
                palette = c("forestgreen"= min(avg.h.data$mean_corr), 
                            "white"= 0, 
                            "firebrick"= max(avg.h.data$mean_corr)),
                show.legend = F) %>% 
    remove_axes()
  
  p <- subplot(lh,rh,sc) %>%
    layout(scene = scene)  %>% 
    layout(title = paste0("predicted activity map for: ", med)) %>%
    config(
      toImageButtonOptions = list(
        format = "svg",
        filename = med,
        width = 1100,
        height = 700
      )
    )
  htmlwidgets::saveWidget(as_widget(p), paste0("data/maps/model-082223/interactive-maps/", med, ".html"))
  gc()
  #################################################################################
}
