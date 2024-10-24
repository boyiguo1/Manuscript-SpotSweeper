library(SpotSweeper)
library(RANN)
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(escheR)
library(patchwork)
library(scran)
library(scater)
library(tidyr)
library(dplyr)


plot_dir <- here("plots","VisiumHD","human_breast")


spe <- readRDS(here("processed-data", "VisiumHD","human_breast", "VisiumHD_HumanBreast_016_bayesspace.rds"))
spe

colnames(colData(spe.sw))
#  [1] "barcode"                       "in_tissue"                    
#  [3] "array_row"                     "array_col"                    
#  [5] "sample_id"                     "sum"                          
#  [7] "detected"                      "subsets_mito_percent"         
#  [9] "subsets_mito_sum"              "sum_discard"                  
# [11] "detected_discard"              "subsets_mito_percent_discard" 
# [13] "sum_log"                       "sum_outliers"                 
# [15] "sum_z"                         "detected_log"                 
# [17] "detected_outliers"             "detected_z"                   
# [19] "subsets_mito_percent_log"      "subsets_mito_percent_outliers"
# [21] "subsets_mito_percent_z"        "BS_k10" 


# ====== Add discards by fixed thresholds ========
spe$detected_threshold <- spe.sw$detected < 500
spe$sum_threshold <- spe.sw$sum < 500
spe$subsets_mito_percent_threshold <- spe.sw$subsets_mito_percent >10

# ======= SpotSweeper =========
spe <- localOutliers(spe, metric = "sum", direction = "lower", log = TRUE)
spe <- localOutliers(spe, metric = "detected", direction = "lower", log = TRUE)
spe <- localOutliers(spe, metric = "subsets_mito_percent", direction = "higher", log = TRUE)

# save
saveRDS(spe, here("processed-data", "VisiumHD", "human_breast", "VisiumHD_HumanBreast_016_bspace_spotsweeper.rds"))

