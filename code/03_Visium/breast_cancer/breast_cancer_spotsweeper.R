library(SpotSweeper)
library(RANN)
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(escheR)
library(patchwork)
library(scran)
library(scater)


plot_dir <- here("plots","Visium", "breast_cancer", "miQC")


spe <- readRDS(here("processed-data", "Visium","breast_cancer","BreastCancer_banksy.rds"))
spe


# ===== Add QC Metrics =====

# get mito genes
rownames(spe) <- rowData(spe)$gene_name
is.mito <- grepl("^MT-", rownames(spe))

# subset to SPE to is mito
spe.mito <- spe[rownames(spe) %in% is.mito,]

# add scuttle qc metrics
spe <- scuttle::addPerCellQC(spe, subsets=list(mito=is.mito))

# ====== Add discards by MADs ========
spe$sum_discard <- isOutlier(spe$sum, nmads=3, type="lower", log=TRUE)
spe$detected_discard <- isOutlier(spe$detected, nmads=3, type="lower", log=TRUE)
spe$subsets_mito_percent_discard <- isOutlier(spe$subsets_mito_percent, nmads=3, type="higher")

# ====== Add discards by fixed thresholds ========
spe$detected_threshold <- spe$detected < 500
spe$sum_threshold <- spe$sum < 500
spe$subsets_mito_percent_threshold <- spe$subsets_mito_percent > 30

# ====== Add discards by miQC ======
library(miQC)

model <- mixtureModel(spe)

png(here(plot_dir, "miQC_mixtureModel.png"), width=4, height=4, units="in", res=300)
plotModel(spe, model) + ggtitle("Breast Cancer")
dev.off()

png(here(plot_dir, "miQC_mixtureModelFiltering.png"), width=4, height=4, units="in", res=300)
plotFiltering(spe, model) + ggtitle("Breast Cancer")
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
spe <- localOutliers(spe, metric = "sum", direction = "lower", log = TRUE)
spe <- localOutliers(spe, metric = "detected", direction = "lower", log = TRUE)
spe <- localOutliers(spe, metric = "subsets_mito_percent", direction = "higher", log = TRUE)

saveRDS(spe, here("processed-data", "Visium","breast_cancer","BreastCancer_banksy_spotsweeper.rds"))


sessionInfo()
# other attached packages:
#  [1] flexmix_2.3-19              lattice_0.22-6             
#  [3] miQC_1.12.0                 scater_1.32.1              
#  [5] scran_1.32.0                scuttle_1.14.0             
#  [7] patchwork_1.2.0             escheR_1.4.0               
#  [9] ggplot2_3.5.1               spatialLIBD_1.16.2         
# [11] SpatialExperiment_1.14.0    SingleCellExperiment_1.26.0
# [13] SummarizedExperiment_1.34.0 Biobase_2.64.0             
# [15] GenomicRanges_1.56.1        GenomeInfoDb_1.40.1        
# [17] IRanges_2.38.1              S4Vectors_0.42.1           
# [19] BiocGenerics_0.50.0         MatrixGenerics_1.16.0      
# [21] matrixStats_1.3.0           here_1.0.1                 
# [23] RANN_2.6.1                  SpotSweeper_1.3.1   