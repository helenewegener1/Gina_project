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
# seurat_obj_QC_filtered_list <- readRDS("08_seurat_QC_filtering/out/seurat_obj_QC_filtered_list.rds")
seurat_obj_QC_filtered_list <- readRDS("08_seurat_QC_filtering/out/seurat_obj_QC_filtered_singlets_list.rds")

sample_names <- names(seurat_obj_QC_filtered_list)

# Check that doublets are removed 
# seurat_obj_QC_filtered_list$`HH117-SILP-INF-PC`$scDblFinder.class %>% table()
# seurat_obj_QC_filtered_list$`HH119-SI-PP-GC-AND-PB-AND-TFH-Pool1`$scDblFinder.class %>% table()

# Load Gina annotation file 
broad_annot_file <- read_excel("00_data/Gene_markers_GL_HW_new.xlsx", sheet = "Very broad level")
detailed_annot_file <- read_excel("00_data/Gene_markers_GL_HW_new.xlsx", sheet = "More detailed level")

# Extract cell marker genes as lists
broad_markers <- broad_annot_file %>% as.list()
broad_markers <- map(broad_markers, ~ .x[!is.na(.x)])
names(broad_markers) <- c("T_cell", "B_cell", "DC", "plasmablast_plasma_cell")

detailed_markers <- detailed_annot_file %>% as.list() %>% na.omit()
detailed_markers <- map(detailed_markers, ~ .x[!is.na(.x)])
names(detailed_markers) <- c("TFH_cell", "Naive_B_cell", "Memory_B_cell", "GC_B_cell", 
                             "Activation_markers", "Immunoglobulin_subsets", "Other_markers")

#### Make sure the markers are in the same format as in the seurat object ###### 

# grep("IGHA", rownames(seurat_obj), value = TRUE, ignore.case = TRUE)

# Load sample to get genes that are in the data.  
seurat_obj <- seurat_obj_QC_filtered_list[["HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"]]

source("09_annotation_pre_integration/script/functions.R")

# Update marker format for the broad markers
broad_markers <- update_marker_names(broad_markers, seurat_obj)

# grep("CD21", rownames(seurat_obj), value = TRUE)

# Update marker format for the detailed markers
detailed_markers <- update_marker_names(detailed_markers, seurat_obj)

############################# Get cell cycle score #############################

s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes # MKI67 in here 

################################################################################

################################ Generate Plots ################################

# Prep to save clustered seurat objects
seurat_obj_clustered_list <- rep(0, length(seurat_obj_QC_filtered_list)) %>% as.list()
names(seurat_obj_clustered_list) <- names(seurat_obj_QC_filtered_list)

# Define number of dimensions for clustering and resolution of clusters. 
# n_dim <- 20 # standard flow
# n_dim <- 30 # SCTransform
res <- 0.3

for (method in c("standard", "SCT")){
  
  for (sample_name in sample_names){
    
    # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
    
    # Define specific seurat object 
    seurat_obj <- seurat_obj_QC_filtered_list[[sample_name]]
    
    # N cells 
    n_cells <- ncol(seurat_obj)
    
    # Create directory for plots of specific sample
    out_dir <- glue("09_seurat_QC_clusters/plot_{method}/{sample_name}")
    dir.create(out_dir, showWarnings = FALSE)
    
    # Caluclate cell cycle scores - Already calculated in QC
    # seurat_obj <- CellCycleScoring(seurat_obj,
    #                                s.features = s.genes,
    #                                g2m.features = g2m.genes)
    
    # seurat_obj$S.Score - continuous 
    # seurat_obj$G2M.Score - continuous 
    # seurat_obj$Phase - discrete 
    
    ############################## Seurat workflow ###############################
    if (method == "standard"){
      n_dim <- 20 # standard flow
      seurat_obj <- NormalizeData(seurat_obj)
      seurat_obj <- FindVariableFeatures(seurat_obj)
      seurat_obj <- ScaleData(seurat_obj)
    } else if (method == "SCT"){
      n_dim <- 30 # SCTransform
      seurat_obj <- SCTransform(seurat_obj)
    }
    
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
    # groups <- c("scDblFinder.class", "Phase")
    groups <- c("Phase")
    
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
      
      FeaturePlot(seurat_obj, features = broad_markers[[markers]], ncol = 3) + 
        plot_annotation(title = glue("{markers}"),
                        caption = glue("N cells: {n_cells}"))
      
      ggsave(glue("{out_dir}/{sample_name}_broad_{markers}.pdf"), width = 14, height = 12)
      
    }
    
    ############################################################################## 
    
    ##################### FeaturePlot with detailed_markers ######################
    for (markers in names(detailed_markers)){
      
      FeaturePlot(seurat_obj, features = detailed_markers[[markers]], ncol = 3) + 
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
  saveRDS(seurat_obj_clustered_list, glue("09_seurat_QC_clusters/out/seurat_obj_clustered_list_singlets_{method}.rds"))

}

rm(seurat_obj_clustered_list, seurat_obj_QC_filtered_list)


################################################################################ 
################################### Doublets ###################################
################################################################################ 

# Load data
seurat_obj_QC_filtered_doublets_list <- readRDS("08_seurat_QC_filtering/out/seurat_obj_QC_filtered_doublets_list.rds")

for (sample_name in sample_names){
  
  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  n_doublets <- seurat_obj_QC_filtered_doublets_list[[sample_name]] %>% ncol()
  
  print(sample_name)
  print(glue("N doublets: {n_doublets}"))
  print("---------------------------------------------")
  
}

rm(seurat_obj_QC_filtered_doublets_list)

################################################################################ 
##################################### All ######################################
################################################################################ 

# Load data
seurat_obj_QC_filtered_all_list <- readRDS("08_seurat_QC_filtering/out/seurat_obj_QC_filtered_list.rds")

seurat_obj_clustered_all_list <- rep(0, length(seurat_obj_QC_filtered_all_list)) %>% as.list()
names(seurat_obj_clustered_all_list) <- names(seurat_obj_QC_filtered_all_list)

method <- "standard"
res <- 0.3

for (sample_name in sample_names){
  
  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  seurat_obj <- seurat_obj_QC_filtered_all_list[[sample_name]]
  
  # Calculate and print doublet percentage
  n_doublets <- count(seurat_obj$scDblFinder.class == "doublet") 
  n_cells <- seurat_obj %>% ncol()
  percentage_doublet <- round(n_doublets/n_cells * 100, 1)
  
  print(sample_name)
  print(glue("Doublet percentage: {percentage_doublet}"))
  print("---------------------------------------------")
  
  # invisible(
  #   suppressWarnings(
  #     suppressMessages({
  if (method == "standard"){
    n_dim <- 20 # standard flow
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
    seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
    seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  } else if (method == "SCT"){
    n_dim <- 30 # SCTransform
    seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)
  }

  DefaultAssay(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  ElbowPlot(seurat_obj)
  seurat_obj <- FindNeighbors(seurat_obj,  dims = 1:n_dim, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = res, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:n_dim, verbose = FALSE)
  
  seurat_obj_clustered_all_list[[sample_name]] <- seurat_obj
  #     })
  #   )
  # )
  
}

saveRDS(seurat_obj_clustered_all_list, glue("09_seurat_QC_clusters/out/seurat_obj_clustered_all_list_{method}.rds"))


################################################################################ 
################################ DC Annotation #################################
################################################################################ 

# method <- "standard"
# seurat_obj_QC_filtered_list <- readRDS(glue("09_seurat_QC_clusters/out/seurat_obj_clustered_list_{method}.rds"))

# Venla annotate after human atlas projection of v.8 object?

sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
seurat_obj <- seurat_obj_QC_filtered_list[[sample_name]]

seurat_obj$scDblFinder.class %>% table()

DimPlot(seurat_obj, label = TRUE, group.by = "RNA_snn_res.0.5") + NoLegend()

cluster.markers <- FindAllMarkers(seurat_obj, group.by = "RNA_snn_res.0.5", only.pos = TRUE)
broad_markers[["DC"]] <- c(
    "ITGAX", "HLA-DRA", "HLA-DRB1",
    "CLEC9A",    # cDC1
    "XCR1",      # cDC1
    "CLEC10A",   # cDC2
    "CD1C",      # cDC2
    "LILRA4",    # pDC
    "CLEC4C",    # pDC
    "CD86", "CCR7"   # activated/migratory DCs
  )

cluster.markers[broad_markers[["DC"]], ]

# Define DC cluster numbers
dc_clusters <- c(2, 4, 8, 14)

# Create a new logical column
seurat_obj$DC_bool <- ifelse(seurat_obj$seurat_clusters %in% dc_clusters, TRUE, FALSE)
table(seurat_obj$DC_bool)

# Ginas annotation 
# Gina_seurat <- readRDS("00_data/Gina_HH117_PP_broadAnn.rds")
# table(Gina_seurat$Celltype)

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
