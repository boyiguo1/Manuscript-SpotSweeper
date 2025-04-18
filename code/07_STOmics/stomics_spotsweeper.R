library(SpotSweeper)
library(RANN)
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(escheR)
library(patchwork)
library(scran)
library(scater)


plot_dir <- here("plots","STOmics", "miQC")


spe <- readRDS(here("processed-data", "STOmics","STOmics_prenatal.rds"))
spe


# ===== Add QC Metrics =====

counts(spe) <- assay(spe)

# get mito genes
is.mito <- grepl("^mt-", rownames(spe))

# subset to SPE to is mito
spe.mito <- spe[rownames(spe) %in% is.mito,]

# add scuttle qc metrics
spe <- scuttle::addPerCellQC(spe, subsets=list(mito=is.mito))

# ====== Add discards by MADs ========
spe$sum_discard <- isOutlier(spe$sum, nmads=2, type="lower", log=TRUE)
spe$detected_discard <- isOutlier(spe$detected, nmads=2, type="lower", log=TRUE)
spe$subsets_mito_percent_discard <- isOutlier(spe$subsets_mito_percent, nmads=2, type="higher")

# ====== Add discards by fixed thresholds ========
spe$detected_threshold <- spe$detected < 1000
spe$sum_threshold <- spe$sum < 1000
spe$subsets_mito_percent_threshold <- spe$subsets_mito_percent > 2.5

# ====== Add discards by miQC ======
library(miQC)

model <- mixtureModel(spe)

png(here(plot_dir, "stomics_miQC_mixtureModel.png"), width=4, height=4, units="in", res=300)
plotModel(spe, model) + ggtitle("STOmics")
dev.off()

png(here(plot_dir, "stomics_miQC_mixtureModelFiltering.png"), width=4, height=4, units="in", res=300)
plotFiltering(spe, model) + ggtitle("STOmics")
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
spe <- localOutliers(spe, n_neighbors = 48, metric = "sum", direction = "lower", log = TRUE)
spe <- localOutliers(spe, n_neighbors = 48, metric = "detected", direction = "lower", log = TRUE)
spe <- localOutliers(spe, n_neighbors = 48, metric = "subsets_mito_percent", direction = "higher", log = TRUE)

saveRDS(spe, here("processed-data", "STOmics","STOmics_prenatal_spotsweeper.rds"))