getwd()

library(SeuratObject)
library(Seurat)
library(dplyr)
library(glue)
library(ggplot2)
library(patchwork)
library(readxl)
library(purrr) # map funciton 

# Load data
seurat_obj_QC_filtered_list <- readRDS("08_seurat_QC_filtering/out/seurat_obj_QC_filtered_list.rds")

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
seurat_obj <- seurat_obj_QC_filtered_list[["HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"]]

source("09_annotation_pre_integration/script/functions.R")

# Update marker format for the broad markers
broad_markers <- update_marker_names(broad_markers, seurat_obj)

# Update marker format for the detailed markers
detailed_markers <- update_marker_names(detailed_markers, seurat_obj)

################################################################################

# Prep to save clustered seurat objects
seurat_obj_clustered_list <- list()

# Make FeaturePlots with the marker genes for each sample
# for (sample_name in names(seurat_obj_QC_filtered_list)){
  
  sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  
  # Define specific seurat object 
  seurat_obj <- seurat_obj_QC_filtered_list[[sample_name]]
  
  # Create directory for plots of specific sample
  out_dir <- glue("09_seurat_QC_clusters/plot/{sample_name}")
  dir.create(out_dir, showWarnings = FALSE)
  
  n_dim <- 30
  res <- 0.3
  n_cells <- ncol(seurat_obj) 
  
  # Seurat workflow
  # seurat_obj <- NormalizeData(seurat_obj)
  # seurat_obj <- FindVariableFeatures(seurat_obj)
  # seurat_obj <- ScaleData(seurat_obj)
  # seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c(""))
  # seurat_obj <- SCTransform(seurat_obj)
  seurat_obj <- SCTransform(seurat_obj)
  DefaultAssay(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)
  ElbowPlot(seurat_obj)
  seurat_obj <- FindNeighbors(seurat_obj,  dims = 1:n_dim)
  seurat_obj <- FindClusters(seurat_obj, resolution = res)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:n_dim)
  
  # Feature plots
  features <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "scDblFinder.score")
  # seurat_obj$scDblFinder.score
  for (feature in features){
    
    FeaturePlot(seurat_obj, features = feature) + 
      labs(
        title = feature, 
        caption = glue("N cells: {n_cells}")
      )
    
    ggsave(glue("{out_dir}/{sample_name}_{feature}.pdf"), width = 8, height = 7)
    
  }
  
  # Doublet Dimplot 
  DimPlot(seurat_obj, group.by = "scDblFinder.class") + 
    labs(
      caption = glue("N cells: {n_cells}")
    )
  ggsave(glue("{out_dir}/{sample_name}_scDblFinder.class.pdf"), width = 8, height = 7)
  
  DimPlot(seurat_obj, label = TRUE, group.by = "seurat_clusters") + NoLegend() + 
    labs(
      title = "Seurat clusters", 
      caption = glue("N cells: {n_cells}\nN dim: {n_dim}\nResolution: {res}")
    )
  ggsave(glue("{out_dir}/{sample_name}_DimPlot.pdf"), width = 8, height = 8)
  
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
  
  # Save clustered seurat object 
  seurat_obj_clustered_list[[sample_name]] <- seurat_obj
  
# }

################################ DC Annotation #################################

sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
seurat_obj <- seurat_obj_clustered_list[[sample_name]]

DimPlot(seurat_obj, label = TRUE) + NoLegend()

cluster.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)
cluster.markers[broad_markers[["DC"]], ]

# Define DC cluster numbers
dc_clusters <- c(7)

# Create a new logical column
seurat_obj$DC_bool <- ifelse(seurat_obj$seurat_clusters %in% dc_clusters, TRUE, FALSE)
table(seurat_obj$DC_bool)

# Ginas annotation 
Gina_seurat <- readRDS("00_data/Gina_HH117_PP_broadAnn.rds")
table(Gina_seurat$Celltype)

# Plot
# DimPlot(seurat_obj, group.by = "seurat_clusters", label = TRUE) + NoLegend() + 
#   DimPlot(seurat_obj, group.by = "DC_bool", label = TRUE, cols = c("grey", "orange")) + NoLegend() 


# Split object in DCs and non-DCs
# First, DCs
seurat_obj_DC_list <- list()

seurat_obj_DC <- subset(seurat_obj, subset = DC_bool == TRUE)
seurat_obj_DC[["ADT"]] <- NULL

seurat_obj_DC_list[[sample_name]] <- seurat_obj_DC

# Then, other
seurat_obj_nonDC_list <- list()

seurat_obj_nonDC <- subset(seurat_obj, subset = DC_bool == FALSE)

seurat_obj_nonDC_list[[sample_name]] <- seurat_obj_nonDC

# Save lists of DC and nonDC seurat objects 
saveRDS(seurat_obj_DC_list, "09_annotation_pre_integration/out/seurat_obj_DC_list.rds")
saveRDS(seurat_obj_nonDC_list, "09_annotation_pre_integration/out/seurat_obj_nonDC_list.rds")


################################## Annotation ##################################
