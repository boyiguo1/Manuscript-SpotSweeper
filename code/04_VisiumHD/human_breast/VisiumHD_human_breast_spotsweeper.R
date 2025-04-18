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


plot_dir <- here("plots","VisiumHD","human_breast", "016")


spe <- readRDS(here("processed-data", "VisiumHD","human_breast", "VisiumHD_HumanBreast_016_bayesspace.rds"))
spe

colnames(colData(spe))
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
spe$detected_threshold <- spe$detected < 500
spe$sum_threshold <- spe$sum < 500
spe$subsets_mito_percent_threshold <- spe$subsets_mito_percent >10



# ====== Add discards by miQC ======
library(miQC)

model <- mixtureModel(spe)

png(here(plot_dir, "miQC_mixtureModel.png"), width=4, height=4, units="in", res=300)
plotModel(spe, model) + ggtitle("HD Breast Cancer")
dev.off()

png(here(plot_dir, "miQC_mixtureModelFiltering.png"), width=4, height=4, units="in", res=300)
plotFiltering(spe, model) + ggtitle("HD Breast Cancer")
dev.off()


# NOTE: miQC does not offer a way to simply flag spots, so I'm importing part of the filter function here to do so manually
library(flexmix)

posterior_cutoff = 0.75
keep_all_below_boundary = TRUE
enforce_left_cutoff = TRUE

metrics <- as.data.frame(colData(spe))

    if (is.null(model)) {
        warning("call 'mixtureModel' explicitly to get stable model features")
        model <- mixtureModel(spe)
    }

    intercept1 <- parameters(model, component = 1)[1]
    intercept2 <- parameters(model, component = 2)[1]
    if (intercept1 > intercept2) {
        compromised_dist <- 1
        intact_dist <- 2
    } else {
        intact_dist <- 1
        compromised_dist <- 2
    }

    post <- posterior(model)
    metrics$prob_compromised <- post[, compromised_dist]
    spe$prob_compromised <- metrics$prob_compromised
    metrics$keep <- metrics$prob_compromised <= posterior_cutoff

    if (sum(metrics$keep)==nrow(metrics)){
        stop("all cells passed posterior probability filtering. One 
              cause of this is the model selecting two near-identical
              distributions. Try rerunning mixtureModel() and/or 
              setting a different random seed.")
    }
    
    if (keep_all_below_boundary == TRUE) {
        predictions <- fitted(model)[, intact_dist]
        metrics$intact_prediction <- predictions
        metrics[metrics$subsets_mito_percent <
                    metrics$intact_prediction, ]$keep <- TRUE
    }

    if (enforce_left_cutoff == TRUE) {
        min_discard <- min(metrics[!metrics$keep, ]$subsets_mito_percent)
        min_index <- which(metrics$subsets_mito_percent == min_discard)[1]
        lib_complexity <- metrics[min_index, ]$detected
        metrics[metrics$detected <= lib_complexity &
                    metrics$subsets_mito_percent >= min_discard, ]$keep <- FALSE
    }

# add "outliers" to SPE
spe$miQC_keep <- !metrics$keep


# ======= SpotSweeper =========
spe <- localOutliers(spe,  n_neighbors = 48, metric = "sum", direction = "lower", log = TRUE)
spe <- localOutliers(spe,  n_neighbors = 48, metric = "detected", direction = "lower", log = TRUE)
spe <- localOutliers(spe,  n_neighbors = 48, metric = "subsets_mito_percent", direction = "higher", log = TRUE)

# save
saveRDS(spe, here("processed-data", "VisiumHD", "human_breast", "VisiumHD_HumanBreast_016_spotsweeper.rds"))

