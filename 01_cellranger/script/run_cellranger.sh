#!/bin/sh
#PBS -W group_list=dtu_00062 -A dtu_00062 # CIIR
#PBS -l nodes=1:ppn=40:thinnode
#PBS -l gpus=0
#PBS -l mem=80GB
#PBS -l walltime=10:00:00
#PBS -e /home/projects/dtu_00062/people/helweg/projects/Gina_project/01_cellranger/script/run_cellranger.err
#PBS -o /home/projects/dtu_00062/people/helweg/projects/Gina_project/01_cellranger/script/run_cellranger.log
#PBS -N run_cellranger

WD="/home/people/helweg/ciir/people/helweg/projects/Gina_project/01_cellranger/out"

module load tools
module load cellranger/9.0.1

# Navigate to your desired output location
cd "$WD"

# Reference from Fred:
# --transcriptome /services/tools/cellranger/refdata/refdata-cellranger-mm10-2.1.0

cellranger count \
  --id=LP_CRAM1 \
  --transcriptome=/home/projects/dtu_00062/people/helweg/LPcDC_cellranger/mouse_reference/refdata-gex-GRCm39-2024-A \
  --fastqs=../fastqs \
  --sample=ERR12552061 \
  --localcores=8 \
  --localmem=64 \
  --create-bam=false

cellranger count \
  --id=LP_CRAM2 \
  --transcriptome=/home/projects/dtu_00062/people/helweg/LPcDC_cellranger/mouse_reference/refdata-gex-GRCm39-2024-A \
  --fastqs=../fastqs \
  --sample=ERR12552062 \
  --localcores=8 \
  --localmem=64 \
  --create-bam=false
