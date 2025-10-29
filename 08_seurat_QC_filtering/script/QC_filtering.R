# setwd("~/ciir/people/helweg/projects/Gina_project/")

library(SeuratObject)
library(Seurat)
library(dplyr)
library(glue)
library(ggplot2)
library(patchwork)

source("08_seurat_QC_filtering/script/functions.R")

# Load data
seurat_obj_list <- readRDS("07_seurat_QC/out/seurat_obj_QC.rds")

# Initialize filtered list
seurat_obj_QC_filtered_list <- list()

# Check samples
names(seurat_obj_list)

# HH117-SI-MILF-INF-HLADR-AND-CD19 looks different than the others

# Doublets out object 
doublet_nFeature_test <- list()

# Extra plots
# for (sample_name in names(seurat_obj_list)){
#   
#   seurat_obj <- seurat_obj_list[[sample_name]]
#   seurat_obj <- pre_filter_pipeline(seurat_obj) # calculate the percent.mt and percent.hb
#   
#   FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
#   FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#   VlnPlot(seurat_obj, features = "percent.hb", layer = "counts")
#   
# }

################################### Rough QC ################################### 

sample_name <- names(seurat_obj_list)[[1]]
seurat_obj <- seurat_obj_list[[sample_name]]

# Doublets out object 
doublet_nFeature_test[[sample_name]] <- table(seurat_obj[[]]$scDblFinder.class, seurat_obj[[]]$nFeature_RNA > 6000)

# Filter on doublets and calculate QC metrics + plot pre_filter plots
seurat_obj <- pre_filter_pipeline(seurat_obj)

# Filter cells based on QC plots
filtering_expr <- expr(nFeature_RNA > 400 & nFeature_RNA < 6000 & percent.mt < 20)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, sample_name = sample_name,
        version = "filtered", filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj, seurat_obj_filtered, n_cells_raw, n_cells_filtered)

################################################################################

sample_name <- names(seurat_obj_list)[[2]]
seurat_obj <- seurat_obj_list[[sample_name]]

# Doublets out object 
doublet_nFeature_test[[sample_name]] <- table(seurat_obj[[]]$DF.classifications, seurat_obj[[]]$nFeature_RNA > 6000)

# Filter on doublets and calculate QC metrics + plot pre_filter plots
seurat_obj <- pre_filter_pipeline(seurat_obj)

# Filter cells based on QC plots
filtering_expr <- expr(nFeature_RNA > 400 & nFeature_RNA < 6000 & percent.mt < 20)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, sample_name = sample_name,
        version = "filtered", filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj, seurat_obj_filtered, n_cells_raw, n_cells_filtered)

################################################################################

sample_name <- names(seurat_obj_list)[[3]]
seurat_obj <- seurat_obj_list[[sample_name]]

# Doublets out object 
doublet_nFeature_test[[sample_name]] <- table(seurat_obj[[]]$DF.classifications, seurat_obj[[]]$nFeature_RNA > 8000)

# Filter on doublets and calculate QC metrics + plot pre_filter plots
seurat_obj <- pre_filter_pipeline(seurat_obj)

# Filter cells based on QC plots
filtering_expr <- expr(nFeature_RNA > 2000 & nFeature_RNA < 8000 & percent.mt < 20)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, sample_name = sample_name, 
        version = "filtered", filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj, seurat_obj_filtered, n_cells_raw, n_cells_filtered)

################################################################################

sample_name <- names(seurat_obj_list)[[4]]
seurat_obj <- seurat_obj_list[[sample_name]]

# Doublets out object 
doublet_nFeature_test[[sample_name]] <- table(seurat_obj[[]]$DF.classifications, seurat_obj[[]]$nFeature_RNA > 8000)

# Filter on doublets and calculate QC metrics + plot pre_filter plots
seurat_obj <- pre_filter_pipeline(seurat_obj)

# Filter cells based on QC plots
filtering_expr <- expr(nFeature_RNA > 400 & nFeature_RNA < 8000 & percent.mt < 20)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, sample_name = sample_name, 
        version = "filtered", filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj, seurat_obj_filtered, n_cells_raw, n_cells_filtered)

################################################################################

sample_name <- names(seurat_obj_list)[[5]]
seurat_obj <- seurat_obj_list[[sample_name]]

# Doublets out object 
doublet_nFeature_test[[sample_name]] <- table(seurat_obj[[]]$DF.classifications, seurat_obj[[]]$nFeature_RNA > 8000)

# Filter on doublets and calculate QC metrics + plot pre_filter plots
seurat_obj <- pre_filter_pipeline(seurat_obj)

# Filter cells based on QC plots
filtering_expr <- expr(nFeature_RNA > 400 & nFeature_RNA < 8000 & percent.mt < 20)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, sample_name = sample_name, 
        version = "filtered", filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj, seurat_obj_filtered, n_cells_raw, n_cells_filtered)

################################################################################

sample_name <- names(seurat_obj_list)[[6]]
seurat_obj <- seurat_obj_list[[sample_name]]

# Doublets out object 
doublet_nFeature_test[[sample_name]] <- table(seurat_obj[[]]$DF.classifications, seurat_obj[[]]$nFeature_RNA > 6000)

# Filter on doublets and calculate QC metrics + plot pre_filter plots
seurat_obj <- pre_filter_pipeline(seurat_obj)

# Filter cells based on QC plots
filtering_expr <- expr(nFeature_RNA > 400 & nFeature_RNA < 6000 & percent.mt < 20)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, sample_name = sample_name, 
        version = "filtered", filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj, seurat_obj_filtered, n_cells_raw, n_cells_filtered)

################################################################################

sample_name <- names(seurat_obj_list)[[7]]
seurat_obj <- seurat_obj_list[[sample_name]]

# Doublets out object 
doublet_nFeature_test[[sample_name]] <- table(seurat_obj[[]]$DF.classifications, seurat_obj[[]]$nFeature_RNA > 7500)

# Filter on doublets and calculate QC metrics + plot pre_filter plots
seurat_obj <- pre_filter_pipeline(seurat_obj)

# Filter cells based on QC plots
filtering_expr <- expr(nFeature_RNA > 400 & nFeature_RNA < 7500 & percent.mt < 20)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)
n_cells_filtered <- ncol(seurat_obj_filtered)

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, sample_name = sample_name, 
        version = "filtered", filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj, seurat_obj_filtered, n_cells_raw, n_cells_filtered)

################################################################################

sample_name <- names(seurat_obj_list)[[8]]
seurat_obj <- seurat_obj_list[[sample_name]]

# Doublets out object 
doublet_nFeature_test[[sample_name]] <- table(seurat_obj[[]]$DF.classifications, seurat_obj[[]]$nFeature_RNA > 6000)

# Filter on doublets and calculate QC metrics + plot pre_filter plots
seurat_obj <- pre_filter_pipeline(seurat_obj)

# Filter cells based on QC plots
filtering_expr <- expr(nFeature_RNA > 400 & nFeature_RNA < 6000 & percent.mt < 20)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)
n_cells_filtered <- ncol(seurat_obj_filtered)

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, sample_name = sample_name, 
        version = "filtered", filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj, seurat_obj_filtered, n_cells_raw, n_cells_filtered)

################################################################################

sample_name <- names(seurat_obj_list)[[9]]
seurat_obj <- seurat_obj_list[[sample_name]]

# Doublets out object 
doublet_nFeature_test[[sample_name]] <- table(seurat_obj[[]]$DF.classifications, seurat_obj[[]]$nFeature_RNA > 7500)

# Filter on doublets and calculate QC metrics + plot pre_filter plots
seurat_obj <- pre_filter_pipeline(seurat_obj)

# Filter cells based on QC plots
filtering_expr <- expr(nFeature_RNA > 400 & nFeature_RNA < 7500 & percent.mt < 20)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)
n_cells_filtered <- ncol(seurat_obj_filtered)

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, sample_name = sample_name, 
        version = "filtered", filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj, seurat_obj_filtered, n_cells_raw, n_cells_filtered)

################################################################################

sample_name <- names(seurat_obj_list)[[10]]
seurat_obj <- seurat_obj_list[[sample_name]]

# Doublets out object 
doublet_nFeature_test[[sample_name]] <- table(seurat_obj[[]]$DF.classifications, seurat_obj[[]]$nFeature_RNA > 6000)

# Filter on doublets and calculate QC metrics + plot pre_filter plots
seurat_obj <- pre_filter_pipeline(seurat_obj)

# Filter cells based on QC plots
filtering_expr <- expr(nFeature_RNA > 400 & nFeature_RNA < 6000 & percent.mt < 20)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)
n_cells_filtered <- ncol(seurat_obj_filtered)

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, sample_name = sample_name, 
        version = "filtered", filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj, seurat_obj_filtered, n_cells_raw, n_cells_filtered)

################################################################################

sample_name <- names(seurat_obj_list)[[11]]

seurat_obj <- seurat_obj_list[[sample_name]]

# Filter on QC metrics from 07
seurat_obj <- subset(seurat_obj, subset = scDblFinder.class == "singlet")

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
filtering_expr <- expr(nFeature_RNA > 400 & nFeature_RNA < 6000 & percent.mt < 20)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)

n_cells_filtered <- ncol(seurat_obj_filtered)

# Doublets out object 
doublet_nFeature_test[[sample_name]] <- table(seurat_obj[[]]$DF.classifications, seurat_obj[[]]$nFeature_RNA > 6000)

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, sample_name = sample_name, 
        version = "filtered", filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj, seurat_obj_filtered, n_cells_raw, n_cells_filtered)

################################################################################

sample_name <- names(seurat_obj_list)[[12]]

seurat_obj <- seurat_obj_list[[sample_name]]

# Filter on QC metrics from 07
seurat_obj <- subset(seurat_obj, subset = scDblFinder.class == "singlet")

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
filtering_expr <- expr(nFeature_RNA > 400 & nFeature_RNA < 6000 & percent.mt < 20)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)

# Doublets out object 
doublet_nFeature_test[[sample_name]] <- table(seurat_obj[[]]$DF.classifications, seurat_obj[[]]$nFeature_RNA > 6000)

n_cells_filtered <- ncol(seurat_obj_filtered)

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, sample_name = sample_name, 
        version = "filtered", filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj, seurat_obj_filtered, n_cells_raw, n_cells_filtered)

################################################################################

sample_name <- names(seurat_obj_list)[[13]]

seurat_obj <- seurat_obj_list[[sample_name]]

# Filter on QC metrics from 07
seurat_obj <- subset(seurat_obj, subset = scDblFinder.class == "singlet")

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
filtering_expr <- expr(nFeature_RNA > 400 & nFeature_RNA < 8000 & percent.mt < 20)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)

n_cells_filtered <- ncol(seurat_obj_filtered)

# Doublets out object 
doublet_nFeature_test[[sample_name]] <- table(seurat_obj[[]]$DF.classifications, seurat_obj[[]]$nFeature_RNA > 8000)

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, sample_name = sample_name, 
        version = "filtered", filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj, seurat_obj_filtered, n_cells_raw, n_cells_filtered)

########################################## Export list of filtered Seurat objects ##########################################

names(seurat_obj_QC_filtered_list)
saveRDS(seurat_obj_QC_filtered_list, "08_seurat_QC_filtering/out/seurat_obj_QC_filtered_list.rds")

# Doublets out object 
doublet_nFeature_test
saveRDS(doublet_nFeature_test, "08_seurat_QC_filtering/out/doublet_nFeature_test.rds")

