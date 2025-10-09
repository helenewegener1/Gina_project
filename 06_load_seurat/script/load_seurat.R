setwd("~/Documents/projects/project_Bcells/Gina_project/")

library(tidyverse)
library(Seurat)
library(glue)

# Sample
sample <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"

# Define path to cellranger output
OUTS_DIR <- glue("05_run_cellranger/out/res_{sample}/outs")

# Step 1: Import Gene Expression (GEX) Data

# Read GEX counts (the default method for 10x data)
gex.data <- Read10X(data.dir = file.path(OUTS_DIR, "filtered_feature_bc_matrix"))

# Create the base Seurat object
# Note: The GEX data is automatically stored in the 'RNA' assay.
seurat_obj <- CreateSeuratObject(
  counts = gex.data$`Gene Expression`,
  project = "Multi_Analysis",
  min.cells = 3,
  min.features = 200
)

# Step 2: Import Antibody Capture (ADT) Data

# Extract the ADT count matrix from the list
adt_counts <- gex.data$`Antibody Capture`

# Add ADT data as a new assay (e.g., "ADT")
seurat_obj[["ADT"]] <- CreateAssayObject(counts = adt_counts)

# Normalize the ADT data (often using the CLR method)
seurat_obj <- NormalizeData(
  seurat_obj, 
  assay = "ADT", 
  normalization.method = "CLR"
)

# Find variable features in the ADT assay (optional but good practice)
seurat_obj <- FindVariableFeatures(
  seurat_obj, 
  selection.method = "vst", 
  assay = "ADT"
)

# Step 3: Import VDJ Data (TCR and BCR)

# Load and filter the VDJ annotation files
tcr_meta <- read.csv(glue("{OUTS_DIR}/vdj_t/filtered_contig_annotations.csv"))
bcr_meta <- read.csv(glue("{OUTS_DIR}/vdj_b/filtered_contig_annotations.csv"))

# We must ensure VDJ data is filtered down to one productive chain per cell 
# and format it for merging (using Seurat's recommended functions or custom scripts)

# --- Use the Seurat utility functions for VDJ data preparation ---
# This step aggregates the contigs (multiple sequences per cell) into
# summary information (e.g., VDJ clonotype, chains, and productivity).

# Add TCR data
seurat_obj <- AddMetaData(
  seurat_obj, 
  metadata = tcr_meta, 
  col.name = "TCR_meta"
)

# Add BCR data
seurat_obj <- AddMetaData(
  seurat_obj, 
  metadata = bcr_meta, 
  col.name = "BCR_meta"
)

# CRITICAL STEP: Add clonotype info (the actual T/B cell ID)
# This usually requires running a specialized function to simplify the contig data.
# E.g., Use the 'clonotype' and 'raw_clonotype_id' columns from the filtered_contig_annotations.csv
# and link them to the cell barcodes (BCR-CELL_BARCODE-1, TCR-CELL_BARCODE-1). 
# Be mindful of the barcode naming convention differences between GEX and VDJ outputs!






