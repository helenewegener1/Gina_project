
# Define working directory 
WD="/home/people/helweg/ciir/people/helweg/projects/Gina_project/"
OUT_DIR="05_run_cellranger/out/"
CONFIG_DIR="04_prep_config/out"

module load tools
module load cellranger/9.0.1

# Navigate to your desired output location
cd "$WD"

# List of all your unique sample identifiers (e.g., HH117-SI-PP-...)
# Use read -a (read into an array) to populate samples_array
read -r -d '' -a samples_array < <(ls "$SAMPLE_DIR" | awk -F '_' '{print $1}' | awk -F '-' '{sub($1 FS, "", $0); print $0}' | sort | uniq)

# Run cellranger multi for all IDs
for ID in "${samples_array[@]}"; do
    
    echo "--- Starting analysis for sample: ${ID} ---"
    echo "${OUT_DIR}/${ID}_Multi_Analysis"
    echo "${CONFIG_DIR}/multi_config_${ID}"

    # cellranger multi \
    #     --id="${OUT_DIR}/${ID}_Multi_Analysis" \
    #     --csv="${CONFIG_DIR}/multi_config_${ID}" \
    #     --localcores=32 \
    #     --localmem=128
    
done
