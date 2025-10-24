getwd()

library(SeuratObject)
library(Seurat)
library(dplyr)
library(glue)
library(ggplot2)
library(patchwork)
library(readxl)
library(scRepertoire)

# Load data
seurat_obj_list <- readRDS("07_seurat_QC/out/seurat_obj_QC_metrics.rds")
# seurat_obj_list <- readRDS("08_seurat_QC_filtering/out/seurat_obj_QC_filtered_list.rds")

# Investigate non-unique chains across contigs which should be filtered on umi. 

# Subset cells with TCR and BCR respectively. 
# Which samples have TCR and BCR data 
tcr_mask <- lapply(names(seurat_obj_list), function(x) {"tcr_v_gene_contig_1" %in% colnames(seurat_obj_list[[x]]@meta.data)}) %>% unlist()
bcr_mask <- lapply(names(seurat_obj_list), function(x) {"bcr_v_gene_contig_1" %in% colnames(seurat_obj_list[[x]]@meta.data)}) %>% unlist()

# Extract seurat objects with TCR and BCR data 
tcr_seurat_obj_list <- seurat_obj_list[tcr_mask]
bcr_seurat_obj_list <- seurat_obj_list[bcr_mask]

##################################### TCR ###################################### 

tcr_double_trouble <- list() 

for (sample_name in names(tcr_seurat_obj_list)){
  
  tcr_meta <- tcr_seurat_obj_list[[sample_name]]@meta.data
  
  colnames(sample_name)
  
  tcr_double_trouble[[sample_name]] <- tcr_meta %>% 
    select(starts_with("tcr_chain"), starts_with("tcr_umi")) %>% 
    mutate(double_trouble = apply(select(., starts_with("tcr_chain")), 1, function(x) {
      x <- x[!is.na(x)]              # remove NAs for comparison, but don't modify the main data
      any(duplicated(x))             # TRUE if any non-NA contigs are duplicated
    })) %>% 
    filter(double_trouble)
  
}

##################################### BCR ###################################### 

bcr_double_trouble <- list() 

for (sample_name in names(bcr_seurat_obj_list)){
  
  bcr_meta <- bcr_seurat_obj_list[[sample_name]]@meta.data
  
  colnames(sample_name)
  
  bcr_double_trouble[[sample_name]] <- bcr_meta %>% 
    select(starts_with("bcr_chain"), starts_with("bcr_umi")) %>% 
    mutate(double_trouble = apply(select(., starts_with("bcr_chain")), 1, function(x) {
      x <- x[!is.na(x)]              # remove NAs for comparison, but don't modify the main data
      any(duplicated(x))             # TRUE if any non-NA contigs are duplicated
    })) %>% 
    filter(double_trouble)
  
}

################################# scRepertorie ################################# 

sample_name <- names(tcr_seurat_obj_list)[[1]]

contigs <- read.csv(glue("05_run_cellranger/out/res_HH117-SI-MILF-INF-HLADR-AND-CD19/outs/per_sample_outs/res_HH117-SI-MILF-INF-HLADR-AND-CD19/web_summary.html//filtered_contig_annotations.csv"))

contig.list <- createHTOContigList(contigs, 
                                   Seurat.Obj, 
                                   group.by = "HTO_maxID")
