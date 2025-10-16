getwd()

library(SeuratObject)
library(Seurat)
library(dplyr)
library(glue)
library(ggplot2)
library(patchwork)
library(readxl)

# Load data
seurat_obj_roughQC_list <- readRDS("07_seurat_roughQC/out/seurat_obj_roughQC_list.rds")

# Load Gina annotation file 
broad_annot_file <- read_excel("00_data/Gene_markers_GL_HW.xlsx", sheet = "Very broad level")
detailed_annot_file <- read_excel("00_data/Gene_markers_GL_HW.xlsx", sheet = "More detailed level")

# Extract cell marker genes as lists
broad_markers <- broad_annot_file %>% as.list()
broad_markers <- map(broad_markers, ~ .x[!is.na(.x)])
names(broad_markers) <- c("CD4_T_cell", "B_cell", "DC", "plasmablast_plasma_cell", "IGHG", "IGHA")

detailed_markers <- detailed_annot_file %>% as.list() %>% na.omit()
detailed_markers <- map(detailed_markers, ~ .x[!is.na(.x)])
names(detailed_markers) <- c("TFH_cell", "Naive_B_cell", "Memory_B_cell", "GC_B_cell")

################## Make sure the markers are in the same format as in the seurat object ################## 

# grep("IGHA", rownames(seurat_obj), value = TRUE, ignore.case = TRUE)

# Load sample to get genes that are in the data.  
seurat_obj <- seurat_obj_roughQC_list[[1]]

source("09_annotation_pre_integration/script/functions.R")

# Update marker format for the broad markers
broad_markers <- update_marker_names(broad_markers, seurat_obj)

# Update marker format for the detailed markers
detailed_markers <- update_marker_names(detailed_markers, seurat_obj)

# Make FeaturePlots with the marker genes for each sample
for (sample_name in names(seurat_obj_roughQC_list)){
  
  # Define specific seurat object 
  seurat_obj <- seurat_obj_roughQC_list[[sample_name]]
  # seurat_obj <- seurat_obj_roughQC_list$`HH119-COLP-PC`
  
  # Create directory for plots of specific sample
  out_dir <- glue("09_annotation_pre_integration/plot/{sample_name}")
  dir.create(out_dir, showWarnings = FALSE)
  
  # Seurat workflow
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)
  ElbowPlot(seurat_obj)
  seurat_obj <- FindNeighbors(seurat_obj,  dims = 1:20)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.4)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:20)
  
  # Number of cells 
  n_cells <- ncol(seurat_obj) 
  
  # Visualize with UMAP stratified by hashtag/follicle
  if ("ADT_maxID" %in% colnames(seurat_obj[[]])) {
    
    DimPlot(seurat_obj, reduction = "umap", group.by = "ADT_maxID") +
      labs(title = "UMAP") +
      theme(legend.text = element_text(size = 8))
    
    ggsave(glue("{out_dir}/{sample_name}_DimPlot_ADT_maxID.pdf"), width = 7, height = 6)
    
  }
  
  # FeaturePlot with broad_markers 
  for (markers in names(broad_markers)){
    
    FeaturePlot(seurat_obj, features = broad_markers[[markers]], ncol = 2) + 
      plot_annotation(title = glue("{markers}"),
                      caption = glue("N cells: {n_cells}"))
    
    ggsave(glue("{out_dir}/{sample_name}_FeaturePlot_broad_{markers}.pdf"), width = 14, height = 12)
    
  }
  
  # FeaturePlot with detailed_markers 
  for (markers in names(detailed_markers)){
    
    FeaturePlot(seurat_obj, features = detailed_markers[[markers]], ncol = 2) + 
      plot_annotation(title = glue("{markers}"), 
                      caption = glue("N cells: {n_cells}"))
    
    ggsave(glue("{out_dir}/{sample_name}_FeaturePlot_detailed_{markers}.pdf"), width = 14, height = 18)
    
  }
  
}








