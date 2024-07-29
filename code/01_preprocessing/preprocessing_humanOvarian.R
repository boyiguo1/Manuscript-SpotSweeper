###################################
# This script was adopted from:
# https://github.com/kinnaryshah/MeanVarBias/blob/main/R/01_preprocessing/preprocessing_humanLobularBreast.R

# dataset: Visium human ovarian cancer
###################################
library(here)
library(SpatialExperiment)
library(nnSVG)
library(scran)
# BiocManager::install("DropletUtils") if not installed

# ---------
# load data
# ---------

spe <- read10xVisium(here("raw-data/humanOvarian/outs"),
                     type = "sparse",   # use sparse (not HDF5) format
                     data = "filtered",
                     images = "lowres", # specify which image(s) to include
                     load = TRUE)      # specify whether or not to load image(s)

# -------------
# preprocessing
# -------------

#row data formatting
rowData(spe)$gene_id <- rownames(spe)
rowData(spe)$gene_name <- rowData(spe)$symbol
rowData(spe)$symbol <- NULL

## Remove genes without enough data
#spe <- filter_genes(spe, 2, 0.2)
#dim(spe)
# > dim(spe)
# [1] 12022  1935

# normalization and log2 transformation
spe <- logNormCounts(spe)


# -----------
# save object
# -----------

fn <- here("processed-data", "spe_humanOvarian_preprocessed.rds")
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
#   [1] nnSVG_1.6.4                 scater_1.30.1               scran_1.30.2                scuttle_1.12.0
# [5] patchwork_1.2.0             escheR_1.2.0                ggplot2_3.5.1               spatialLIBD_1.14.1
# [9] SpatialExperiment_1.12.0    SingleCellExperiment_1.24.0 SummarizedExperiment_1.32.0 Biobase_2.62.0
# [13] GenomicRanges_1.54.1        GenomeInfoDb_1.38.8         IRanges_2.36.0              S4Vectors_0.40.2
# [17] BiocGenerics_0.48.1         MatrixGenerics_1.14.0       matrixStats_1.3.0           here_1.0.1
# [21] RANN_2.6.1                  SpotSweeper_0.99.4
#
# loaded via a namespace (and not attached):
#   [1] later_1.3.2                   BiocIO_1.12.0                 bitops_1.0-7
# [4] filelock_1.0.3                fields_15.2                   R.oo_1.26.0
# [7] tibble_3.2.1                  XML_3.99-0.16.1               lifecycle_1.0.4
# [10] edgeR_4.0.16                  doParallel_1.0.17             rprojroot_2.0.4
# [13] lattice_0.22-5                MASS_7.3-60.0.1               magrittr_2.0.3
# [16] sass_0.4.9                    limma_3.58.1                  plotly_4.10.4
# [19] jquerylib_0.1.4               yaml_2.3.8                    remotes_2.5.0
# [22] metapod_1.10.1                httpuv_1.6.15                 spam_2.10-0
# [25] sessioninfo_1.2.2             pbapply_1.7-2                 cowplot_1.1.3
# [28] DBI_1.2.2                     RColorBrewer_1.1-3            golem_0.4.1
# [31] maps_3.4.2                    abind_1.4-5                   zlibbioc_1.48.2
# [34] rvest_1.0.4                   R.utils_2.12.3                purrr_1.0.2
# [37] RCurl_1.98-1.14               rappdirs_0.3.3                GenomeInfoDbData_1.2.11
# [40] ggrepel_0.9.5                 irlba_2.3.5.1                 terra_1.7-71
# [43] dqrng_0.3.2                   DelayedMatrixStats_1.24.0     DropletUtils_1.22.0
# [46] codetools_0.2-19              DelayedArray_0.28.0           xml2_1.3.6
# [49] DT_0.33                       tidyselect_1.2.1              farver_2.1.1
# [52] ScaledMatrix_1.10.0           viridis_0.6.5                 shinyWidgets_0.8.6
# [55] BiocFileCache_2.10.2          GenomicAlignments_1.38.2      jsonlite_1.8.8
# [58] BiocNeighbors_1.20.2          iterators_1.0.14              foreach_1.5.2
# [61] tools_4.3.3                   Rcpp_1.0.12                   glue_1.7.0
# [64] gridExtra_2.3                 SparseArray_1.2.4             HDF5Array_1.30.1
# [67] dplyr_1.1.4                   withr_3.0.0                   BiocManager_1.30.22
# [70] fastmap_1.1.1                 rhdf5filters_1.14.1           bluster_1.12.0
# [73] fansi_1.0.6                   digest_0.6.35                 rsvd_1.0.5
# [76] R6_2.5.1                      mime_0.12                     colorspace_2.1-0
# [79] RSQLite_2.3.6                 R.methodsS3_1.8.2             config_0.3.2
# [82] utf8_1.2.4                    tidyr_1.3.1                   generics_0.1.3
# [85] data.table_1.15.4             rtracklayer_1.62.0            httr_1.4.7
# [88] htmlwidgets_1.6.4             S4Arrays_1.2.1                pkgconfig_2.0.3
# [91] gtable_0.3.5                  rdist_0.0.5                   blob_1.2.4
# [94] BRISC_1.0.5                   selectr_0.4-2                 XVector_0.42.0
# [97] htmltools_0.5.8.1             dotCall64_1.1-1               fastqcr_0.1.3
# [100] scales_1.3.0                  png_0.1-8                     attempt_0.3.1
# [103] rstudioapi_0.16.0             rjson_0.2.21                  curl_5.2.1
# [106] rhdf5_2.46.1                  cachem_1.0.8                  stringr_1.5.1
# [109] BiocVersion_3.18.1            parallel_4.3.3                vipor_0.4.7
# [112] AnnotationDbi_1.64.1          restfulr_0.0.15               pillar_1.9.0
# [115] grid_4.3.3                    vctrs_0.6.5                   promises_1.3.0
# [118] BiocSingular_1.18.0           dbplyr_2.5.0                  beachmat_2.18.1
# [121] xtable_1.8-4                  cluster_2.1.6                 beeswarm_0.4.0
# [124] paletteer_1.6.0               magick_2.8.3                  Rsamtools_2.18.0
# [127] cli_3.6.2                     locfit_1.5-9.9                compiler_4.3.3
# [130] rlang_1.1.3                   crayon_1.5.2                  labeling_0.4.3
# [133] rematch2_2.1.2                spatialEco_2.0-2              ggbeeswarm_0.7.2
# [136] stringi_1.8.3                 viridisLite_0.4.2             BiocParallel_1.36.0
# [139] munsell_0.5.1                 Biostrings_2.70.3             lazyeval_0.2.2
# [142] Matrix_1.6-5                  ExperimentHub_2.10.0          benchmarkme_1.0.8
# [145] sparseMatrixStats_1.14.0      bit64_4.0.5                   Rhdf5lib_1.24.2
# [148] KEGGREST_1.42.0               statmod_1.5.0                 shiny_1.8.1.1
# [151] interactiveDisplayBase_1.40.0 AnnotationHub_3.10.1          igraph_2.0.3
# [154] memoise_2.0.1                 bslib_0.7.0                   benchmarkmeData_1.0.4
# [157] bit_4.0.5
