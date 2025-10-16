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

lapply(names(seurat_obj_roughQC_list), "contig_1_v_gene" %in% colnames(seurat_obj_roughQC_list[[x]]@meta.data))
