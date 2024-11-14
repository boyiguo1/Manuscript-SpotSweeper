#!/bin/bash
#SBATCH -p transfer
#SBATCH -c 1
#SBATCH --mem=3G
#SBATCH --job-name=get_Xenium5k
#SBATCH --output=./logs/get_Xenium5k.out
#SBATCH --error=./logs/get_Xenium5k.err
#SBATCH --open-mode=append


# originally copied from /dcs05/lieber/lcolladotor/Visium_HD_DLPFC_pilot_LIBD4100/Visium_HD_DLPFC_pilot/

repo_dir=$(git rev-parse --show-toplevel)
dest_dir=$repo_dir/raw-data/Xenium

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"


# =============================================================================
#   Sample 1: Mouse Brain
# =============================================================================
mkdir -p $dest_dir/MouseBrain_5k/input
mkdir $dest_dir/MouseBrain_5k/output

#   All input files 
cd $dest_dir/MouseBrain_5k/input
curl -O https://cf.10xgenomics.com/samples/xenium/3.0.0/Xenium_Prime_Mouse_Brain_Coronal_FF/Xenium_Prime_Mouse_Brain_Coronal_FF_gene_panel.json

#   All output files
cd $dest_dir/MouseBrain_5k/output
curl -O https://cf.10xgenomics.com/samples/xenium/3.0.0/Xenium_Prime_Mouse_Brain_Coronal_FF/Xenium_Prime_Mouse_Brain_Coronal_FF_outs.zip
curl -O https://cf.10xgenomics.com/samples/xenium/3.0.0/Xenium_Prime_Mouse_Brain_Coronal_FF/Xenium_Prime_Mouse_Brain_Coronal_FF_xe_outs.zip


#   Now extract and unzip any tar.gz files
for tar_file in $(ls *.tar.gz); do
    tar -xzf $tar_file
done


# =============================================================================
#   Sample 2: Human Breast Cancer
# =============================================================================
mkdir -p $dest_dir/HumanBreast_5k/input
mkdir $dest_dir/HumanBreast_5k/output

#  Input
cd $dest_dir/HumanBreast_5k/output
curl -O https://cf.10xgenomics.com/samples/xenium/3.0.0/Xenium_Prime_Breast_Cancer_FFPE/Xenium_Prime_Breast_Cancer_FFPE_gene_panel.json


#  Output
cd $dest_dir/HumanBreast_5k/output
curl -O https://s3-us-west-2.amazonaws.com/10x.files/samples/xenium/3.0.0/Xenium_Prime_Breast_Cancer_FFPE/Xenium_Prime_Breast_Cancer_FFPE_outs.zip
curl -O https://cf.10xgenomics.com/samples/xenium/3.0.0/Xenium_Prime_Breast_Cancer_FFPE/Xenium_Prime_Breast_Cancer_FFPE_xe_outs.zip


#   Now extract and unzip any tar.gz files
for tar_file in $(ls *.tar.gz); do
    tar -xzf $tar_file
done


echo "**** Job ends ****"
date