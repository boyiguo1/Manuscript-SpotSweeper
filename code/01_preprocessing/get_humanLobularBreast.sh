#!/bin/bash
#SBATCH --job-name=get_humanOvarian
#SBATCH --output=./logs/get_humanLobularBreast.out
#SBATCH --error=./logs/get_humanLobularBreast.err
#SBATCH --mem-per-cpu=5G --cpus-per-task=1
#SBATCH --array=1


# make directory for the data first time
# mkdir outputs/raw
cd ..
cd ..
cd raw-data

mkdir humanLobularBreast
cd humanLobularBreast

mkdir outs
cd outs

curl -O https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_BreastCancer/Parent_Visium_Human_BreastCancer_raw_feature_bc_matrix.tar.gz
curl -O https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_BreastCancer/Parent_Visium_Human_BreastCancer_spatial.tar.gz

gunzip Parent_Visium_Human_BreastCancer_raw_feature_bc_matrix.tar.gz
tar -xvf Parent_Visium_Human_BreastCancer_raw_feature_bc_matrix.tar

gunzip Parent_Visium_Human_BreastCancer_spatial.tar.gz
tar -xvf Parent_Visium_Human_BreastCancer_spatial.tar

# clean up
rm Parent_Visium_Human_BreastCancer_raw_feature_bc_matrix.tar
rm Parent_Visium_Human_BreastCancer_spatial.tar
