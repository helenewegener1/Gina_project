#!/bin/sh
#PBS -W group_list=dtu_00062 -A dtu_00062 # CIIR
#PBS -l nodes=1:ppn=32:thinnode
#PBS -l gpus=0
#PBS -l mem=32GB
#PBS -l walltime=08:00:00
#PBS -e /home/projects/dtu_00062/people/helweg/projects/Gina_project/04_prep_config/script/prep_config.err
#PBS -o /home/projects/dtu_00062/people/helweg/projects/Gina_project/04_prep_config/script/prep_config.log
#PBS -N prep_config

# Define working directory 
WD="/home/people/helweg/ciir/people/helweg/projects/Gina_project/04_prep_config"
SAMPLE_DIR="/home/projects/dtu_00062/data/KU09/FASTQ_ku09_mkfastq/outs/fastq_path/HKL3YDSXF"

# Navigate to your desired output location
cd "$WD"

# List of all your unique sample identifiers (e.g., HH117-SI-PP-...)
# Use read -a (read into an array) to populate samples_array
read -r -d '' -a samples_array < <(ls "$SAMPLE_DIR" | awk -F '_' '{print $1}' | awk -F '-' '{sub($1 FS, "", $0); print $0}' | sort | uniq)


for ID in "${SAMPLE_IDS[@]}"; do
    # 1. Generate the final CSV file for this sample ID
    TEMP_CSV="${WD}/${CONFIG_DIR}/multi_config_GEX_ADT_BCR_TCR_${ID}.csv"
    
    # Use sed to replace the placeholder "SAMPLE_PREFIX" with the actual ID
    sed "s/SAMPLE_PREFIX/${ID}/g" multi_config_GEX_ADT_BCR_TCR_template.csv > "${TEMP_CSV}"
    
    
    echo "--- Finished config file for sample: ${ID} ---"
    
done


