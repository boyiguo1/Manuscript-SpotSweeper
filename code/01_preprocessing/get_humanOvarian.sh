#!/bin/bash
#SBATCH --job-name=get_humanOvarian
#SBATCH --output=./logs/get_humanOvarian.out
#SBATCH --error=./logs/get_humanOvarian.err
#SBATCH --mem-per-cpu=5G --cpus-per-task=1
#SBATCH --array=1

# List of sample IDs
samples=("GSM6506110" "GSM6506111" "GSM6506112" "GSM6506113" "GSM6506114" "GSM6506115" "GSM6506116" "GSM6506117")

# Base directory for raw data
cd ..
cd ..
base_dir=$(pwd)
data_dir="${base_dir}/raw-data/humanOvarian/outs"

# Create base directories if they don't exist
mkdir -p "${data_dir}"


# Loop through each sample ID and download files
for i in "${!samples[@]}"; do
    sample="${samples[$i]}"
    sp_num=$((i + 1))  # SP number (SP1, SP2, ..., SP8)
    sp_label="SP${sp_num}"
    sample_dir="${data_dir}/${sample}"
    mkdir -p "${sample_dir}/filtered_feature_bc_matrix"

    # Download spatial data
    curl -o "${sample_dir}/spatial.zip" "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6506nnn/${sample}/suppl/${sample}_${sp_label}_spatial.zip"
    unzip "${sample_dir}/spatial.zip" -d "${sample_dir}"

    # Download filtered feature bc matrix data
    curl -o "${sample_dir}/filtered_feature_bc_matrix/barcodes.tsv.gz" "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6506nnn/${sample}/suppl/${sample}_${sp_label}_barcodes.tsv.gz"
    curl -o "${sample_dir}/filtered_feature_bc_matrix/features.tsv.gz" "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6506nnn/${sample}/suppl/${sample}_${sp_label}_features.tsv.gz"
    curl -o "${sample_dir}/filtered_feature_bc_matrix/matrix.mtx.gz" "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6506nnn/${sample}/suppl/${sample}_${sp_label}_matrix.mtx.gz"

    # Clean up
    rm "${sample_dir}/spatial.zip"
done
