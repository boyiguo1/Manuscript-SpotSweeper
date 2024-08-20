#!/bin/bash

#SBATCH -p transfer
#SBATCH -c 1
#SBATCH --mem=3G
#SBATCH -t 1-00:00
#SBATCH --job-name=01_download_10x_examples
#SBATCH -o ../../processed-data/01_first_look/01_download_10x_examples.log
#SBATCH -e ../../processed-data/01_first_look/01_download_10x_examples.log
#SBATCH --open-mode=append

# originally copied from /dcs05/lieber/lcolladotor/Visium_HD_DLPFC_pilot_LIBD4100/Visium_HD_DLPFC_pilot/

repo_dir=$(git rev-parse --show-toplevel)
dest_dir=$repo_dir/processed-data/01_first_look

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

################################################################################
#   Sample 1: Visium_HD_Mouse_Brain
################################################################################

mkdir -p $dest_dir/Visium_HD_Mouse_Brain/input
mkdir $dest_dir/Visium_HD_Mouse_Brain/output

#   All input files except FASTQs
cd $dest_dir/Visium_HD_Mouse_Brain/input
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Brain/Visium_HD_Mouse_Brain_image.tif
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Brain/Visium_HD_Mouse_Brain_tissue_image.tif
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Brain/Visium_HD_Mouse_Brain_alignment_file.json
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Brain/Visium_HD_Mouse_Brain_probe_set.csv

#   All output files
cd $dest_dir/Visium_HD_Mouse_Brain/output
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Brain/Visium_HD_Mouse_Brain_web_summary.html
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Brain/Visium_HD_Mouse_Brain_cloupe_008um.cloupe
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Brain/Visium_HD_Mouse_Brain_feature_slice.h5
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Brain/Visium_HD_Mouse_Brain_metrics_summary.csv
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Brain/Visium_HD_Mouse_Brain_molecule_info.h5
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Brain/Visium_HD_Mouse_Brain_spatial.tar.gz
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Brain/Visium_HD_Mouse_Brain_square_002um_outputs.tar.gz
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Brain/Visium_HD_Mouse_Brain_square_008um_outputs.tar.gz
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Brain/Visium_HD_Mouse_Brain_square_016um_outputs.tar.gz

#   Now extract and unzip any tar.gz files
for tar_file in $(ls *.tar.gz); do
    tar -xzf $tar_file
done

# ################################################################################
# #   Sample 2: Visium_HD_Human_Lung_Cancer
# ################################################################################

# mkdir -p $dest_dir/Visium_HD_Human_Lung_Cancer/input
# mkdir $dest_dir/Visium_HD_Human_Lung_Cancer/output

#   All input files except FASTQs
cd $dest_dir/Visium_HD_Human_Lung_Cancer/input
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Lung_Cancer/Visium_HD_Human_Lung_Cancer_image.tif
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Lung_Cancer/Visium_HD_Human_Lung_Cancer_tissue_image.tif
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Lung_Cancer/Visium_HD_Human_Lung_Cancer_probe_set.csv

#   All output files
cd $dest_dir/Visium_HD_Human_Lung_Cancer/output
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Lung_Cancer/Visium_HD_Human_Lung_Cancer_web_summary.html
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Lung_Cancer/Visium_HD_Human_Lung_Cancer_cloupe_008um.cloupe
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Lung_Cancer/Visium_HD_Human_Lung_Cancer_feature_slice.h5
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Lung_Cancer/Visium_HD_Human_Lung_Cancer_metrics_summary.csv
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Lung_Cancer/Visium_HD_Human_Lung_Cancer_molecule_info.h5
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Lung_Cancer/Visium_HD_Human_Lung_Cancer_spatial.tar.gz
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Lung_Cancer/Visium_HD_Human_Lung_Cancer_square_002um_outputs.tar.gz
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Lung_Cancer/Visium_HD_Human_Lung_Cancer_square_008um_outputs.tar.gz
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Lung_Cancer/Visium_HD_Human_Lung_Cancer_square_016um_outputs.tar.gz

#   Now extract and unzip any tar.gz files
for tar_file in $(ls *.tar.gz); do
    tar -xzf $tar_file
done

# ################################################################################
# #   Sample 3: Visium_HD_Mouse_Small_Intestine
# ################################################################################

# mkdir -p $dest_dir/Visium_HD_Mouse_Small_Intestine/input
# mkdir $dest_dir/Visium_HD_Mouse_Small_Intestine/output
#
# #   All input files except FASTQs
# cd $dest_dir/Visium_HD_Mouse_Small_Intestine/input
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Small_Intestine/Visium_HD_Mouse_Small_Intestine_image.tif
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Small_Intestine/Visium_HD_Mouse_Small_Intestine_tissue_image.btf
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Small_Intestine/Visium_HD_Mouse_Small_Intestine_probe_set.csv
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Small_Intestine/Visium_HD_Mouse_Small_Intestine_slide_file.vlf
#
# #   All output files
# cd $dest_dir/Visium_HD_Mouse_Small_Intestine/output
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Small_Intestine/Visium_HD_Mouse_Small_Intestine_web_summary.html
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Small_Intestine/Visium_HD_Mouse_Small_Intestine_cloupe_008um.cloupe
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Small_Intestine/Visium_HD_Mouse_Small_Intestine_feature_slice.h5
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Small_Intestine/Visium_HD_Mouse_Small_Intestine_metrics_summary.csv
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Small_Intestine/Visium_HD_Mouse_Small_Intestine_molecule_info.h5
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Small_Intestine/Visium_HD_Mouse_Small_Intestine_spatial.tar.gz
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Small_Intestine/Visium_HD_Mouse_Small_Intestine_square_002um_binned_outputs.tar.gz
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Small_Intestine/Visium_HD_Mouse_Small_Intestine_square_008um_binned_outputs.tar.gz
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Small_Intestine/Visium_HD_Mouse_Small_Intestine_square_016um_binned_outputs.tar.gz
#
# #   Now extract and unzip any tar.gz files
# for tar_file in $(ls *.tar.gz); do
#     tar -xzf $tar_file
# done

################################################################################
#   Sample 4: Visium_HD_Human_Colon_Cancer
################################################################################
#
# mkdir -p $dest_dir/Visium_HD_Human_Colon_Cancer/input
# mkdir $dest_dir/Visium_HD_Human_Colon_Cancer/output
#
# echo "Downloading input files for Visium_HD_Human_Colon_Cancer..."
#
# #   All input files except FASTQs
# cd $dest_dir/Visium_HD_Human_Colon_Cancer/input
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_image.tif
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_tissue_image.btf
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_alignment_file.json
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_probe_set.csv
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_slide_file.vlf
#
# echo "Downloading output files for Visium_HD_Human_Colon_Cancer..."
#
# #   All output files
# cd $dest_dir/Visium_HD_Human_Colon_Cancer/output
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_web_summary.html
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_cloupe_008um.cloupe
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_feature_slice.h5
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_metrics_summary.csv
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_molecule_info.h5
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_spatial.tar.gz
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_square_002um_outputs.tar.gz
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_square_008um_outputs.tar.gz
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_square_016um_outputs.tar.gz
#
# #   Now extract and unzip any tar.gz files
# for tar_file in $(ls *.tar.gz); do
# tar -xzf $tar_file
# done

echo "**** Job ends ****"
date
