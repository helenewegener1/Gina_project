getwd()

library(SeuratObject)
library(Seurat)
library(dplyr)
library(glue)
library(ggplot2)
library(patchwork)

# Load data
seurat_obj_roughQC_list <- readRDS("07_seurat_roughQC/out/seurat_obj_roughQC_list.rds")

sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"

seurat_obj <- seurat_obj_roughQC_list[[sample_name]]

# Seurat workflow
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)

ElbowPlot(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj,  dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.4)
seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:20)

# Visualize with UMAP stratified by hashtag/follicle
if ("ADT_maxID" %in% colnames(seurat_obj[[]])) {
  
  DimPlot(seurat_obj, reduction = "umap", group.by = "ADT_maxID") +
    labs(title = "UMAP") +
    theme(legend.text = element_text(size = 8))
  
  ggsave(glue("09_annotation_pre_integration/plot/{sample_name}_DimPlot_ADT_maxID.pdf"), width = 7, height = 6)
  
}

FeaturePlot(seurat_obj, features = c("CD19")) + 
  plot_annotation(title = "B-cells")
ggsave(glue("09_annotation_pre_integration/plot/{sample_name}_FeaturePlot_CD19.pdf"), width = 7, height = 6)

FeaturePlot(seurat_obj, features = c("CD19", "HLA-DRB1")) + 
  plot_annotation(title = "GC: Germinal Center (B-cells)")
# grep("HLA-DRB", rownames(seurat_obj), value = TRUE, ignore.case = TRUE)
ggsave(glue("09_annotation_pre_integration/plot/{sample_name}_FeaturePlot_CG.pdf"), width = 12, height = 6)

FeaturePlot(seurat_obj, features = c("CXCR5")) + 
  plot_annotation(title = "TFH: T Follicular Helper (T-cells)")
ggsave(glue("09_annotation_pre_integration/plot/{sample_name}_FeaturePlot_TFH.pdf"), width = 7, height = 6)

FeaturePlot(seurat_obj, features = c("CD19", "CD38", "PRDM1"), ncol = 2) + 
  plot_annotation(title = "PB: Plasmablast")
ggsave(glue("09_annotation_pre_integration/plot/{sample_name}_FeaturePlot_PB.pdf"), width = 12, height = 12)

FeaturePlot(seurat_obj, features = c("PRDM1", "IRF4"), ncol = 2) + 
  plot_annotation(title = "PC: Plasma Cell")
ggsave(glue("09_annotation_pre_integration/plot/{sample_name}_FeaturePlot_PC.pdf"), width = 12, height = 6)

FeaturePlot(seurat_obj, features = c("HLA-DRA", "HLA-DRB5", "HLA-DRB6", "HLA-DRB1"), ncol = 2) + 
  plot_annotation(title = "HLADR: Human Leukocyte Antigen – (MHC Class II molecule)")
ggsave(glue("09_annotation_pre_integration/plot/{sample_name}_FeaturePlot_HLADR.pdf"), width = 14, height = 12)

# grep("HLA-DR", rownames(seurat_obj), value = TRUE, ignore.case = TRUE)

  

