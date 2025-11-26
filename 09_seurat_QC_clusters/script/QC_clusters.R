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

#### Make sure the markers are in the same format as in the seurat object ###### 

# grep("IGHA", rownames(seurat_obj), value = TRUE, ignore.case = TRUE)

# Load sample to get genes that are in the data.  
seurat_obj <- seurat_obj_QC_filtered_list[["HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"]]

source("09_annotation_pre_integration/script/functions.R")

# Update marker format for the broad markers
broad_markers <- update_marker_names(broad_markers, seurat_obj)

# Update marker format for the detailed markers
detailed_markers <- update_marker_names(detailed_markers, seurat_obj)

############################# Get cell cycle score #############################

s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes

################################################################################

################################ Generate Plots ################################

# Prep to save clustered seurat objects
seurat_obj_clustered_list <- list()

# Define number of dimensions for clustering and resolution of clusters. 
# n_dim <- 20 # standard flow
n_dim <- 30 # SCTransform
res <- 0.3

for (sample_name in names(seurat_obj_QC_filtered_list)){
  
  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  
  # Define specific seurat object 
  seurat_obj <- seurat_obj_QC_filtered_list[[sample_name]]
  
  # N cells 
  n_cells <- ncol(seurat_obj)
  
  # Create directory for plots of specific sample
  out_dir <- glue("09_seurat_QC_clusters/plot/{sample_name}")
  dir.create(out_dir, showWarnings = FALSE)
  
  # Caluclate cell cycle scores
  seurat_obj <- CellCycleScoring(seurat_obj,
                                 s.features = s.genes,
                                 g2m.features = g2m.genes)
  
  # seurat_obj$S.Score - continuous 
  # seurat_obj$G2M.Score - continuous 
  # seurat_obj$Phase - discrete 
  
  ############################## Seurat workflow ###############################
  # seurat_obj <- NormalizeData(seurat_obj)
  # seurat_obj <- FindVariableFeatures(seurat_obj)
  # seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- SCTransform(seurat_obj)
  DefaultAssay(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)
  ElbowPlot(seurat_obj)
  seurat_obj <- FindNeighbors(seurat_obj,  dims = 1:n_dim)
  seurat_obj <- FindClusters(seurat_obj, resolution = res)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:n_dim)
  
  ################## Feature plots with continuous features ####################
  features <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", 
                "scDblFinder.score", "S.Score", "G2M.Score")
  for (feature in features){
    
    FeaturePlot(seurat_obj, features = feature) + 
      labs(
        title = feature, 
        caption = glue("N cells: {n_cells}")
      )
    
    ggsave(glue("{out_dir}/{sample_name}_{feature}.pdf"), width = 8, height = 7)
    
  }
  
  ############################################################################## 
  
  ###################### Dimplots with discrete features #######################
  groups <- c("scDblFinder.class", "Phase")
  
  for (group in groups){
    
    DimPlot(seurat_obj, group.by = group) + 
      labs(
        title = group, 
        caption = glue("N cells: {n_cells}")
      )
    
    ggsave(glue("{out_dir}/{sample_name}_{group}.pdf"), width = 8, height = 7)
    
  }
  
  ############################################################################## 
  
  ####################### FeaturePlot with broad_markers ####################### 
  for (markers in names(broad_markers)){
    
    FeaturePlot(seurat_obj, features = broad_markers[[markers]], ncol = 2) + 
      plot_annotation(title = glue("{markers}"),
                      caption = glue("N cells: {n_cells}"))
    
    ggsave(glue("{out_dir}/{sample_name}_broad_{markers}.pdf"), width = 14, height = 12)
    
  }
  
  ############################################################################## 
  
  ##################### FeaturePlot with detailed_markers ######################
  for (markers in names(detailed_markers)){
    
    FeaturePlot(seurat_obj, features = detailed_markers[[markers]], ncol = 2) + 
      plot_annotation(title = glue("{markers}"), 
                      caption = glue("N cells: {n_cells}"))
    
    ggsave(glue("{out_dir}/{sample_name}_detailed_{markers}.pdf"), width = 14, height = 18)
    
  }
  
  ############################################################################## 
  
  ################################# Clusters ################################# 
  DimPlot(seurat_obj, label = TRUE, group.by = "seurat_clusters") + NoLegend() + 
    labs(
      title = "Seurat clusters", 
      caption = glue("N cells: {n_cells}\nN dim: {n_dim}\nResolution: {res}")
    )
  ggsave(glue("{out_dir}/{sample_name}_clusters.pdf"), width = 8, height = 8)
  
  ############################################################################## 
  
  # Save clustered seurat object 
  seurat_obj_clustered_list[[sample_name]] <- seurat_obj
  
}

# Save object 
saveRDS(seurat_obj_clustered_list, "09_seurat_QC_clusters/out/seurat_obj_clustered_list.rds")

################################################################################ 
################################### Doublets ###################################

# Functions for plotting
source("09_seurat_QC_clusters/script/functions.R")

# Add columns with co-expression of markers from different cell types (very likely doublets)
# Define marker sets
B_markers <- c("MS4A1","CD79A","CD79B","CD19") 
T_markers <- c("CD3D","CD3E","CD3G","TRAC") 
Myeloid_markers <- c("LYZ","S100A8","S100A9","CTSS","FCGR3A") 
Plasma_markers <- c("SDC1","MZB1","XBP1","PRDM1")

marker_pairs <- combn(c("percent_B", "percent_T", "percent_Myeloid", "percent_Plasma"), 2, simplify = FALSE)

for (sample_name in names(seurat_obj_clustered_list)){
  
  sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  seurat_obj <- seurat_obj_clustered_list[[sample_name]]
  
  table(seurat_obj$scDblFinder.class)
  
  # Create directory for plots of specific sample
  out_dir <- glue("09_seurat_QC_clusters/plot/{sample_name}/doublet")
  dir.create(out_dir, showWarnings = FALSE)
  
  doublet_N_genes(seurat_obj, sample_name)
  
  seurat_obj$percent_B       <- PercentageFeatureSet(seurat_obj, features = B_markers)
  seurat_obj$percent_T       <- PercentageFeatureSet(seurat_obj, features = T_markers)
  seurat_obj$percent_Myeloid <- PercentageFeatureSet(seurat_obj, features = Myeloid_markers)
  seurat_obj$percent_Plasma  <- PercentageFeatureSet(seurat_obj, features = Plasma_markers)
  
  for (pair in marker_pairs) {
    
    marker_1 <- pair[1]   # "B"
    marker_2 <- pair[2]   # "T"
    
    # marker_1 <- "percent_B"
    # marker_2 <- "percent_T"     
    
    doublet_dual_lineages(seurat_obj, sample_name, marker_1 = marker_1, marker_2 = marker_2)

  }
  
}



################################################################################ 
################################ DC Annotation #################################

# Venla annotate after human atlas projection of v.8 object?

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
