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

# get folder names in raw-data
samples <- list.dirs(here("raw-data", "humanOvarian"), full.names = TRUE, recursive = FALSE)
sample_id <- list.dirs(here("raw-data", "humanOvarian"), full.names = FALSE, recursive = FALSE)



spe <- read10xVisium(samples=samples,
                     sample_id=sample_id,
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
# R version 4.4.0 Patched (2024-05-22 r86590)
# Platform: x86_64-conda-linux-gnu
# Running under: Rocky Linux 9.4 (Blue Onyx)

# Matrix products: default
# BLAS:   /jhpce/shared/community/core/conda_R/4.4/R/lib64/R/lib/libRblas.so 
# LAPACK: /jhpce/shared/community/core/conda_R/4.4/R/lib64/R/lib/libRlapack.so;  LAPACK version 3.12.0

# locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

# time zone: US/Eastern
# tzcode source: system (glibc)

# attached base packages:
# [1] stats4    stats     graphics  grDevices datasets  utils     methods  
# [8] base     

# other attached packages:
#  [1] scran_1.32.0                scuttle_1.14.0             
#  [3] nnSVG_1.8.0                 SpatialExperiment_1.14.0   
#  [5] SingleCellExperiment_1.26.0 SummarizedExperiment_1.34.0
#  [7] Biobase_2.64.0              GenomicRanges_1.56.1       
#  [9] GenomeInfoDb_1.40.1         IRanges_2.38.1             
# [11] S4Vectors_0.42.1            BiocGenerics_0.50.0        
# [13] MatrixGenerics_1.16.0       matrixStats_1.3.0          
# [15] here_1.0.1                 

# loaded via a namespace (and not attached):
#  [1] rjson_0.2.21              rhdf5_2.48.0             
#  [3] lattice_0.22-6            rhdf5filters_1.16.0      
#  [5] tools_4.4.0               parallel_4.4.0           
#  [7] cluster_2.1.6             R.oo_1.26.0              
#  [9] pkgconfig_2.0.3           BiocNeighbors_1.22.0     
# [11] Matrix_1.7-0              sparseMatrixStats_1.16.0 
# [13] dqrng_0.4.1               lifecycle_1.0.4          
# [15] GenomeInfoDbData_1.2.12   compiler_4.4.0           
# [17] statmod_1.5.0             bluster_1.14.0           
# [19] codetools_0.2-20          crayon_1.5.3             
# [21] R.utils_2.12.3            BiocParallel_1.38.0      
# [23] DelayedArray_0.30.1       limma_3.60.4             
# [25] BRISC_1.0.5               magick_2.8.4             
# [27] abind_1.4-5               metapod_1.12.0           
# [29] locfit_1.5-9.10           rsvd_1.0.5               
# [31] BiocSingular_1.20.0       rprojroot_2.0.4          
# [33] grid_4.4.0                cli_3.6.3                
# [35] SparseArray_1.4.8         magrittr_2.0.3           
# [37] S4Arrays_1.4.1            edgeR_4.2.1              
# [39] DelayedMatrixStats_1.26.0 UCSC.utils_1.0.0         
# [41] XVector_0.44.0            httr_1.4.7               
# [43] DropletUtils_1.24.0       igraph_2.0.3             
# [45] RANN_2.6.1                R.methodsS3_1.8.2        
# [47] ScaledMatrix_1.12.0       beachmat_2.20.0          
# [49] HDF5Array_1.32.1          pbapply_1.7-2            
# [51] irlba_2.3.5.1             rlang_1.1.4              
# [53] Rcpp_1.0.13               rdist_0.0.5              
# [55] jsonlite_1.8.8            Rhdf5lib_1.26.0          
# [57] R6_2.5.1                  zlibbioc_1.50.0  