
# Define working directory 
WD="/home/people/helweg/ciir/people/helweg/projects/Gina_project"
SAMPLE_DIR="/home/projects/dtu_00062/data/KU09/FASTQ_ku09_mkfastq/outs/fastq_path/HKL3YDSXF"

# Navigate to your desired output location
cd "$WD"

# Each unique samples (without GET, ADT, TCR and BCR)
# samples=$(ls $SAMPLE_DIR | awk -F '_' '{print $1}' | awk -F '-' '{sub($1 FS, "", $0); print $0}' | sort | uniq)

# Use read -a (read into an array) to populate samples_array
read -r -d '' -a samples_array < <(ls "$SAMPLE_DIR" | awk -F '_' '{print $1}' | awk -F '-' '{sub($1 FS, "", $0); print $0}' | sort | uniq)


for ID in "${samples_array[@]}"; do

    echo $ID

    
    echo "--- Finished rearrangement for sample: ${ID} ---"
    
done

