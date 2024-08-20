import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import os
import bin2cell as b2c

root = os.getcwd()

processed_path = "../../processed-data/VisiumHD"
source_image_path = "../../raw-data/VisiumHD/Visium_HD_Mouse_Brain/input/Visium_HD_Mouse_Brain_tissue_image.tif"

um2_path = "../../raw-data/VisiumHD/Visium_HD_Mouse_Brain/output/binned_outputs/square_002um/"
um8_path = "../../raw-data/VisiumHD/Visium_HD_Mouse_Brain/output/binned_outputs/square_008um/"
um16_path = "../../raw-data/VisiumHD/Visium_HD_Mouse_Brain/output/binned_outputs/square_016um/"

# read in the data and make names unique
adata_2um = b2c.read_visium(um2_path, source_image_path = source_image_path)
adata_2um.write(processed_path + "/anndats_2um_bins.h5ad")

adata_8um = b2c.read_visium(um8_path, source_image_path = source_image_path)
adata_8um.write(processed_path + "/anndats_8um_bins.h5ad")

adata_16um = b2c.read_visium(um16_path, source_image_path = source_image_path)
adata_16um.write(processed_path + "/anndats_16um_bins.h5ad")