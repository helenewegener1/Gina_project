#!/bin/sh
#PBS -W group_list=dtu_00062 -A dtu_00062 # CIIR
#PBS -l nodes=1:ppn=32:thinnode
#PBS -l gpus=0
#PBS -l mem=128GB
#PBS -l walltime=24:00:00
#PBS -e /home/projects/dtu_00062/people/helweg/projects/Gina_project/02_cellranger/script/run_cellranger.err
#PBS -o /home/projects/dtu_00062/people/helweg/projects/Gina_project/02_cellranger/script/run_cellranger.log
#PBS -N run_cellranger

# Define working directory 
WD="/home/people/helweg/ciir/people/helweg/projects/Gina_project"

module load tools
module load cellranger/9.0.1

# Navigate to your desired output location
cd "$WD"

# List of all your unique sample identifiers (e.g., HH117-SI-PP-...)
SAMPLE_IDS=("HH117-SI-PP-nonINF..." "HH119-SI-PP-GC-AND-PB...")

for ID in "${SAMPLE_IDS[@]}"; do
    # 1. Generate the final CSV file for this sample ID
    TEMP_CSV="${ID}_config.csv"
    
    # Use sed to replace the placeholder "SAMPLE_PREFIX" with the actual ID
    sed "s/SAMPLE_PREFIX/${ID}/g" config_template.csv > "${TEMP_CSV}"
    
    echo "--- Starting analysis for sample: ${ID} ---"
    
    # 2. Run cellranger multi with the dynamically generated CSV
    cellranger multi \
        --id="../out/${ID}_Multi_Analysis" \
        --csv="${TEMP_CSV}" \
        --localcores=32 \
        --localmem=128G
    
    echo "--- Finished analysis for sample: ${ID} ---"
done

# Reference from Fred:
# --transcriptome /services/tools/cellranger/refdata/refdata-cellranger-mm10-2.1.0

# cellranger count \
#   --id=LP_CRAM1 \
#   --transcriptome=/home/projects/dtu_00062/people/helweg/LPcDC_cellranger/mouse_reference/refdata-gex-GRCm39-2024-A \
#   --fastqs=../fastqs \
#   --sample=ERR12552061 \
#   --localcores=8 \
#   --localmem=64 \
#   --create-bam=false
