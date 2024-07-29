###################################
# This script was adapted from:
# https://github.com/kinnaryshah/MeanVarBias/blob/main/R/01_preprocessing/preprocessing_humanLobularBreast.R
#
# dataset: Visium human breast cancer
# https://www.10xgenomics.com/datasets/human-breast-cancer-whole-transcriptome-analysis-1-standard-1-2-0
###################################
library(scran)
library(SpatialExperiment)
library(here)
# BiocManager::install("DropletUtils") if not installed


# ---------
# load data
# ---------

spe <-  read10xVisium(here("raw-data/humanLobularBreast/outs"),
                      type = "sparse",   # use sparse (not HDF5) format
                      data = "raw",
                      images = "lowres", # specify which image(s) to include
                      load = TRUE)      # specify whether or not to load image(s)

# -------------
# preprocessing
# -------------

#row data formatting
rowData(spe)$gene_id <- rownames(spe)
rowData(spe)$gene_name <- rowData(spe)$symbol
rowData(spe)$symbol <- NULL

#remove 667 spots not in tissue
spe <- spe[, colData(spe)$in_tissue == 1]

## remove genes without enough data
#spe <- nnSVG::filter_genes(spe, 2, 0.2)
#dim(spe)
# [1] 12624  4325

# create logcounts matrix
spe <- logNormCounts(spe)


# -----------
# save object
# -----------

fn <- here("processed-data", "spe_humanLobularBreast_preprocessed.rds")
saveRDS(spe, file = fn)

# -----------
# session information
# -----------

sessionInfo()
# R version 4.3.3 (2024-02-29)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Sonoma 14.5
#
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
#
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#
# time zone: America/New_York
# tzcode source: internal
#
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#   [1] scater_1.30.1               scran_1.30.2                scuttle_1.12.0              patchwork_1.2.0
# [5] escheR_1.2.0                ggplot2_3.5.1               spatialLIBD_1.14.1          SpatialExperiment_1.12.0
# [9] SingleCellExperiment_1.24.0 SummarizedExperiment_1.32.0 Biobase_2.62.0              GenomicRanges_1.54.1
# [13] GenomeInfoDb_1.38.8         IRanges_2.36.0              S4Vectors_0.40.2            BiocGenerics_0.48.1
# [17] MatrixGenerics_1.14.0       matrixStats_1.3.0           here_1.0.1                  RANN_2.6.1
# [21] SpotSweeper_0.99.4
#
# loaded via a namespace (and not attached):
#   [1] later_1.3.2                   BiocIO_1.12.0                 bitops_1.0-7
# [4] filelock_1.0.3                fields_15.2                   tibble_3.2.1
# [7] XML_3.99-0.16.1               lifecycle_1.0.4               edgeR_4.0.16
# [10] doParallel_1.0.17             rprojroot_2.0.4               lattice_0.22-5
# [13] MASS_7.3-60.0.1               magrittr_2.0.3                sass_0.4.9
# [16] limma_3.58.1                  plotly_4.10.4                 jquerylib_0.1.4
# [19] yaml_2.3.8                    remotes_2.5.0                 metapod_1.10.1
# [22] httpuv_1.6.15                 spam_2.10-0                   sessioninfo_1.2.2
# [25] pbapply_1.7-2                 cowplot_1.1.3                 DBI_1.2.2
# [28] RColorBrewer_1.1-3            golem_0.4.1                   maps_3.4.2
# [31] abind_1.4-5                   zlibbioc_1.48.2               purrr_1.0.2
# [34] RCurl_1.98-1.14               rappdirs_0.3.3                GenomeInfoDbData_1.2.11
# [37] ggrepel_0.9.5                 irlba_2.3.5.1                 terra_1.7-71
# [40] dqrng_0.3.2                   DelayedMatrixStats_1.24.0     codetools_0.2-19
# [43] DelayedArray_0.28.0           DT_0.33                       tidyselect_1.2.1
# [46] farver_2.1.1                  ScaledMatrix_1.10.0           viridis_0.6.5
# [49] shinyWidgets_0.8.6            BiocFileCache_2.10.2          GenomicAlignments_1.38.2
# [52] jsonlite_1.8.8                BiocNeighbors_1.20.2          iterators_1.0.14
# [55] foreach_1.5.2                 tools_4.3.3                   Rcpp_1.0.12
# [58] glue_1.7.0                    gridExtra_2.3                 SparseArray_1.2.4
# [61] dplyr_1.1.4                   withr_3.0.0                   BiocManager_1.30.22
# [64] fastmap_1.1.1                 bluster_1.12.0                fansi_1.0.6
# [67] digest_0.6.35                 rsvd_1.0.5                    R6_2.5.1
# [70] mime_0.12                     colorspace_2.1-0              RSQLite_2.3.6
# [73] config_0.3.2                  utf8_1.2.4                    tidyr_1.3.1
# [76] generics_0.1.3                data.table_1.15.4             rtracklayer_1.62.0
# [79] httr_1.4.7                    htmlwidgets_1.6.4             S4Arrays_1.2.1
# [82] pkgconfig_2.0.3               gtable_0.3.5                  rdist_0.0.5
# [85] blob_1.2.4                    BRISC_1.0.5                   XVector_0.42.0
# [88] htmltools_0.5.8.1             nnSVG_1.6.4                   dotCall64_1.1-1
# [91] scales_1.3.0                  png_0.1-8                     attempt_0.3.1
# [94] rstudioapi_0.16.0             rjson_0.2.21                  curl_5.2.1
# [97] cachem_1.0.8                  BiocVersion_3.18.1            parallel_4.3.3
# [100] vipor_0.4.7                   AnnotationDbi_1.64.1          restfulr_0.0.15
# [103] pillar_1.9.0                  grid_4.3.3                    vctrs_0.6.5
# [106] promises_1.3.0                BiocSingular_1.18.0           dbplyr_2.5.0
# [109] beachmat_2.18.1               xtable_1.8-4                  cluster_2.1.6
# [112] beeswarm_0.4.0                paletteer_1.6.0               magick_2.8.3
# [115] Rsamtools_2.18.0              cli_3.6.2                     locfit_1.5-9.9
# [118] compiler_4.3.3                rlang_1.1.3                   crayon_1.5.2
# [121] labeling_0.4.3                rematch2_2.1.2                spatialEco_2.0-2
# [124] ggbeeswarm_0.7.2              viridisLite_0.4.2             BiocParallel_1.36.0
# [127] munsell_0.5.1                 Biostrings_2.70.3             lazyeval_0.2.2
# [130] Matrix_1.6-5                  ExperimentHub_2.10.0          benchmarkme_1.0.8
# [133] sparseMatrixStats_1.14.0      bit64_4.0.5                   KEGGREST_1.42.0
# [136] statmod_1.5.0                 shiny_1.8.1.1                 interactiveDisplayBase_1.40.0
# [139] AnnotationHub_3.10.1          igraph_2.0.3                  memoise_2.0.1
# [142] bslib_0.7.0                   benchmarkmeData_1.0.4         bit_4.0.5
