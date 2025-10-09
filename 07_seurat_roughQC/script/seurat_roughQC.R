setwd("~/ciir/people/helweg/projects/Gina_project/")

library(SeuratObject)
library(Seurat)
library(dplyr)
library(glue)
library(ggplot2)
library(patchwork)

source("07_seurat_roughQC/script/functions.R")

# Load data
seurat_obj_list <- readRDS("06_load_seurat/out/seurat_obj_list.rds")

# Initialize filtered list
seurat_obj_roughQC_list <- list()

################################### Rough QC ################################### 

sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"

seurat_obj <- seurat_obj_list[[sample_name]]

# GINA: WHAT DO YOU WANT TO FILTER ON?
# Calculate QC metrics
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RPS|^RPL")
seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^HBA|^HBB")

n_cells_raw <- ncol(seurat_obj) 

# Plot QC metrics in violin plots
plot_qc(seurat_obj = seurat_obj, 
        sample_name = sample_name, 
        n_cells = n_cells_raw, 
        version = "raw", 
        filtering = "")

# Extra plots 
FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
VlnPlot(seurat_obj, features = "percent.hb", layer = "counts")

# Filter cells based on QC plots
filtering_expr <- expr(nFeature_RNA > 1500 & nFeature_RNA < 6000 & percent.mt < 6)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)

n_cells_filtered <- ncol(seurat_obj_filtered)

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, 
        sample_name = sample_name, 
        n_cells = n_cells_filtered, 
        version = "filtered", 
        filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_roughQC_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj, seurat_obj_filtered, n_cells_raw, n_cells_filtered)


################################################################################

sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"


########################################## Export list of filtered Seurat objects ##########################################

names(seurat_obj_roughQC_list)
saveRDS(seurat_obj_roughQC_list, "02_roughQC/out/seurat_obj_roughQC_list.rds")



