
#library(remotes)
#remotes::install_github("pachterlab/SpatialFeatureExperiment", ref = "devel")

library(here)
library(SpatialFeatureExperiment)
library(SpatialExperiment)
library(Voyager)

# define directories
raw_dir <- here("raw-data", "VisiumHD", "Visium_HD_mouse_brain")
processed_dir <- here("processed-data", "VisiumHD")
plot_dir <- here("plots","VisiumHD", "mouse_brain")


hd_dir <-  here("raw-data", "VisiumHD", "Visium_HD_Mouse_Brain", "output", "binned_outputs")

# load VisiumHD dat
sfe.008 <- read10xVisiumSFE(samples=list.files(hd_dir),
                  dirs = hd_dir,
                  bin_size = c("8"), # this defines which of 1:3 resolutions to load
                  type = "HDF5", # Note, "sparse" -> takes longer to load
                  data = "filtered", # spots under tissue
                  images = c("lowres"), # for now low res. image only
                  add_Graph = FALSE # Note, if VisiumHD this can take time for 2 or 8µm res.
                  )

sfe.016 <- read10xVisiumSFE(samples=list.files(hd_dir),
                  dirs = hd_dir,
                  bin_size = c("16"), # this defines which of 1:3 resolutions to load
                  type = "HDF5", # Note, "sparse" -> takes longer to load
                  data = "filtered", # spots under tissue
                  images = c("lowres"), # for now low res. image only
                  add_Graph = FALSE # Note, if VisiumHD this can take time for 2 or 8µm res.
                  )

sfe.002 <- read10xVisiumSFE(samples=list.files(hd_dir),
                  dirs = hd_dir,
                  bin_size = c("2"), # this defines which of 1:3 resolutions to load
                  type = "HDF5", # Note, "sparse" -> takes longer to load
                  data = "filtered", # spots under tissue
                  images = c("lowres"), # for now low res. image only
                  add_Graph = FALSE # Note, if VisiumHD this can take time for 2 or 8µm res.
                  )


# ======== Convert to SpatialExperiment ========

spe <- function(sfe){
  sfe$xcoord <- sfe$array_col
  sfe$ycoord <- sfe$array_row
  
  spe <- SpatialExperiment(assay=counts(sfe),
                           colData=colData(sfe),
                           spatialCoordsNames=c("xcoord", "ycoord"))
  
  rowData(spe) <- rowData(sfe)
  counts(spe) <- assay(sfe)
  
  return(spe)
}

spe.008 <- spe(sfe.008)
spe.016 <- spe(sfe.016)
spe.002 <- spe(sfe.002)
spe.002
# class: SpatialExperiment 
# dim: 19059 6296688 
# metadata(0):
# assays(2): '' counts
# rownames(19059): ENSMUSG00000051951 ENSMUSG00000025900 ...
#   ENSMUSG00000064370 ENSMUSG00000051412
# rowData names(1): symbol
# colnames(6296688): s_002um_02448_01644-1 s_002um_00700_02130-1 ...
#   s_002um_00292_01011-1 s_002um_00865_01560-1
# colData names(5): barcode in_tissue array_row array_col sample_id
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : xcoord ycoord
# imgData names(0):

# save
saveRDS(spe.008, file=here("processed-data", "VisiumHD", "VisiumHD_MouseBrain_008.rds"))
saveRDS(spe.016, file=here("processed-data", "VisiumHD", "VisiumHD_MouseBrain_016.rds"))
saveRDS(spe.002, file=here("processed-data", "VisiumHD", "VisiumHD_MouseBrain_002.rds"))