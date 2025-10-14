setwd("~/ciir/people/helweg/projects/Gina_project/")

library(tidyverse)
library(Seurat)
library(glue)

# Samples
samples <- list.files("05_run_cellranger/out/") %>% str_split_i("_", 2)

# Prep out list
seurat_obj_list <- list()

for (sample in samples){
  
  # sample <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  # sample <- samples[[1]]
  
  # Define path to cellranger output
  # The files in the per_sample_outs directory have been demultiplexed to single samples.
  # Read more about output of cellrange multi: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-3p-outputs-cellplex
  OUTS_DIR <- glue("05_run_cellranger/out/res_{sample}/outs/per_sample_outs/res_{sample}")
  
  # Read GEX counts (the default method for 10x data)
  gex.data <- Read10X(data.dir = glue("{OUTS_DIR}/count/sample_filtered_feature_bc_matrix"))
  
  # Create the base Seurat object
  
  ######################## If antibody data available ######################## 
  if (length(gex.data) == 2 & "Gene Expression" %in% names(gex.data)){ 
    
    
    # Import Gene Expression (GEX) Data 
    seurat_obj <- CreateSeuratObject(
      counts = gex.data$`Gene Expression`,
      project = "Multi_Analysis", 
      min.cells = 3,
      min.features = 200
    )
    
    # Import Antibody Capture (ADT) Data 
    
    # Extract the ADT count matrix from the list
    adt_counts <- gex.data$`Antibody Capture`
    
    # Get the list of cell barcodes that exist in the GEX-based Seurat object
    gex_cells <- colnames(seurat_obj)
    
    # Subset the ADT matrix to include ONLY those GEX cells
    adt_counts_aligned <- adt_counts[, gex_cells]
    
    # Add ADT data as a new assay (e.g., "ADT")
    seurat_obj[["ADT"]] <- CreateAssayObject(counts = adt_counts_aligned) 
    
    # Normalize the ADT data (often using the CLR method)
    seurat_obj <- NormalizeData(seurat_obj, assay = "ADT", normalization.method = "CLR")
    
    # Find variable features in the ADT assay (optional but good practice)
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", assay = "ADT")
    
    
    ################### If only gene expression data is available ################### 
  } else { 
    
    # Import Gene Expression (GEX) Data 
    seurat_obj <- CreateSeuratObject(
      counts = gex.data,
      project = "Multi_Analysis", 
      min.cells = 3,
      min.features = 200
    )
    
  }
  
  #################### Import VDJ Annotation Files (TCR and BCR) ##################### 
  
 
  
  # Check if BCR data is available 
  if ("vdj_t" %in% list.files(OUTS_DIR)){
    tcr_meta <- read.csv(glue("{OUTS_DIR}/vdj_t/filtered_contig_annotations.csv"))
  } 
  
  if ("vdj_b" %in% list.files(OUTS_DIR)){
    bcr_meta <- read.csv(glue("{OUTS_DIR}/vdj_b/filtered_contig_annotations.csv"))
  } 

  # Investigate vdj data
  tcr_meta$barcode %>% length() / 2
  tcr_meta$barcode %>% unique() %>% length()

  bcr_meta$barcode %>% length() / 2
  bcr_meta$barcode %>% unique() %>% length()
  
  bcr_meta$contig_id %>% length()
  bcr_meta$contig_id %>% unique() %>% length()
  
  seurat_obj[[]] %>% nrow()
  
  cells <- rownames(seurat_obj[[]])
  t_cells <- tcr_meta$barcode %>% unique()
  b_cells <- bcr_meta$barcode %>% unique()
  
  table(tcr_meta$is_cell, useNA = "ifany")
  table(bcr_meta$is_cell, useNA = "ifany")
  
  table(t_cells %in% b_cells)
  table(t_cells %in% cells)
  table(b_cells %in% cells)
  b_cells[!b_cells %in% cells]
  
  
  tcr_meta %>% pivot_wider(names_from = contig_id, values_from = v_gene) %>% 
    colnames()
  
  # Test wrangling
  
  # Can we safely remove barcode from contig ids
  table(tcr_meta$barcode == str_split_i(tcr_meta$contig_id, "_", 1))
  table(bcr_meta$barcode == str_split_i(bcr_meta$contig_id, "_", 1))
  
  # GINA: WHICH COLUMNS ARE WE INTERSTED IN?
  tcr_meta_clean <- tcr_meta %>% select(barcode, contig_id, chain, v_gene, d_gene, j_gene, c_gene) %>% 
    mutate(contig_id = str_remove_all(contig_id, glue("{barcode}_"))) %>% 
    pivot_wider(names_from = contig_id, 
                values_from = c(chain, v_gene, d_gene, j_gene, c_gene)) %>% view()

  

  # # We must ensure VDJ data is filtered down to one productive chain per cell 
  # # and format it for merging (using Seurat's recommended functions or custom scripts)
  # 
  # # --- Use the Seurat utility functions for VDJ data preparation ---
  # # This step aggregates the contigs (multiple sequences per cell) into
  # # summary information (e.g., VDJ clonotype, chains, and productivity).
  # 
  # # Add TCR data
  # seurat_obj <- AddMetaData(
  #   seurat_obj, 
  #   metadata = tcr_meta, 
  #   col.name = "TCR_meta"
  # )
  # 
  # # Add BCR data
  # seurat_obj <- AddMetaData(
  #   seurat_obj, 
  #   metadata = bcr_meta, 
  #   col.name = "BCR_meta"
  # )
  # 
  # CRITICAL STEP: Add clonotype info (the actual T/B cell ID)
  # This usually requires running a specialized function to simplify the contig data.
  # E.g., Use the 'clonotype' and 'raw_clonotype_id' columns from the filtered_contig_annotations.csv
  # and link them to the cell barcodes (BCR-CELL_BARCODE-1, TCR-CELL_BARCODE-1). 
  # Be mindful of the barcode naming convention differences between GEX and VDJ outputs!
  
  # Add seurat object to list 
  seurat_obj_list[[sample]] <- seurat_obj
  
}

######################## Export list of seurat objects ######################### 

saveRDS(seurat_obj_list, "06_load_seurat/out/seurat_obj_list.rds")




