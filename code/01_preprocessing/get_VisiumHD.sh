#!/bin/bash
#SBATCH -p transfer
#SBATCH -c 1
#SBATCH --mem=3G
#SBATCH --job-name=get_VisiumHD
#SBATCH --output=./logs/get_VisiumHD.out
#SBATCH --error=./logs/get_VisiumHD.err
#SBATCH --open-mode=append


# originally copied from /dcs05/lieber/lcolladotor/Visium_HD_DLPFC_pilot_LIBD4100/Visium_HD_DLPFC_pilot/

repo_dir=$(git rev-parse --show-toplevel)
dest_dir=$repo_dir/raw-data/VisiumHD

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

# mkdir -p $dest_dir/Visium_HD_Mouse_Brain/input
# mkdir $dest_dir/Visium_HD_Mouse_Brain/output

# #   All input files except FASTQs
# cd $dest_dir/Visium_HD_Mouse_Brain/input
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Brain/Visium_HD_Mouse_Brain_image.tif
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Brain/Visium_HD_Mouse_Brain_tissue_image.tif
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Brain/Visium_HD_Mouse_Brain_alignment_file.json
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Brain/Visium_HD_Mouse_Brain_probe_set.csv

# #   All output files
# cd $dest_dir/Visium_HD_Mouse_Brain/output
# #curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Brain/Visium_HD_Mouse_Brain_web_summary.html
# #curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Brain/Visium_HD_Mouse_Brain_cloupe_008um.cloupe
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Brain/Visium_HD_Mouse_Brain_feature_slice.h5
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Brain/Visium_HD_Mouse_Brain_metrics_summary.csv
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Brain/Visium_HD_Mouse_Brain_molecule_info.h5
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Brain/Visium_HD_Mouse_Brain_spatial.tar.gz
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Brain/Visium_HD_Mouse_Brain_binned_outputs.tar.gz

# #   Now extract and unzip any tar.gz files
# for tar_file in $(ls *.tar.gz); do
#     tar -xzf $tar_file
# done

# ################################################################################
# #   Sample 2: Visium_HD_Human_Lung_Cancer
# ################################################################################

# mkdir -p $dest_dir/Visium_HD_Human_Lung_Cancer/input
# mkdir $dest_dir/Visium_HD_Human_Lung_Cancer/output

# #   All input files except FASTQs
# cd $dest_dir/Visium_HD_Human_Lung_Cancer/input
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Lung_Cancer/Visium_HD_Human_Lung_Cancer_image.tif
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Lung_Cancer/Visium_HD_Human_Lung_Cancer_tissue_image.tif
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Lung_Cancer/Visium_HD_Human_Lung_Cancer_probe_set.csv

# #   All output files
# cd $dest_dir/Visium_HD_Human_Lung_Cancer/output
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Lung_Cancer/Visium_HD_Human_Lung_Cancer_web_summary.html
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Lung_Cancer/Visium_HD_Human_Lung_Cancer_cloupe_008um.cloupe
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Lung_Cancer/Visium_HD_Human_Lung_Cancer_feature_slice.h5
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Lung_Cancer/Visium_HD_Human_Lung_Cancer_metrics_summary.csv
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Lung_Cancer/Visium_HD_Human_Lung_Cancer_molecule_info.h5
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Lung_Cancer/Visium_HD_Human_Lung_Cancer_spatial.tar.gz
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Lung_Cancer/Visium_HD_Human_Lung_Cancer_binned_outputs.tar.gz

# #   Now extract and unzip any tar.gz files
# for tar_file in $(ls *.tar.gz); do
#     tar -xzf $tar_file
# done

################################################################################
#   Sample 3: Visium_HD_Human_Breast_Cancer
################################################################################

mkdir -p $dest_dir/Visium_HD_Human_Breast/input
mkdir $dest_dir/Visium_HD_Human_Breast/output

#   All input files except FASTQs
cd $dest_dir/Visium_HD_Human_Breast/input
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Human_Breast_Cancer_Fresh_Frozen/Visium_HD_Human_Breast_Cancer_Fresh_Frozen_image.tif
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Human_Breast_Cancer_Fresh_Frozen/Visium_HD_Human_Breast_Cancer_Fresh_Frozen_tissue_image.tif
# curl -O https://s3-us-west-2.amazonaws.com/10x.files/samples/spatial-exp/3.1.1/Visium_HD_Human_Breast_Cancer_Fresh_Frozen/Visium_HD_Human_Breast_Cancer_Fresh_Frozen_fastqs.tar
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Human_Breast_Cancer_Fresh_Frozen/Visium_HD_Human_Breast_Cancer_Fresh_Frozen_probe_set.csv

#   All output files
cd $dest_dir/Visium_HD_Human_Breast/output
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Human_Breast_Cancer_Fresh_Frozen/Visium_HD_Human_Breast_Cancer_Fresh_Frozen_web_summary.html
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Human_Breast_Cancer_Fresh_Frozen/Visium_HD_Human_Breast_Cancer_Fresh_Frozen_cloupe_008um.cloupe
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Human_Breast_Cancer_Fresh_Frozen/Visium_HD_Human_Breast_Cancer_Fresh_Frozen_feature_slice.h5
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Human_Breast_Cancer_Fresh_Frozen/Visium_HD_Human_Breast_Cancer_Fresh_Frozen_metrics_summary.csv
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Human_Breast_Cancer_Fresh_Frozen/Visium_HD_Human_Breast_Cancer_Fresh_Frozen_molecule_info.h5
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Human_Breast_Cancer_Fresh_Frozen/Visium_HD_Human_Breast_Cancer_Fresh_Frozen_spatial.tar.gz
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Human_Breast_Cancer_Fresh_Frozen/Visium_HD_Human_Breast_Cancer_Fresh_Frozen_binned_outputs.tar.gz


#   Now extract and unzip any tar.gz files
for tar_file in $(ls *.tar.gz); do
    tar -xzf $tar_file
done

echo "**** Job ends ****"
date