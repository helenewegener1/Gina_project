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

sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
seurat_obj <- seurat_obj_list[[sample_name]]

t_contigs <- read.csv(glue("05_run_cellranger/out/res_{sample_name}/outs/per_sample_outs/res_{sample_name}/vdj_t/filtered_contig_annotations.csv"))
b_contigs <- read.csv(glue("05_run_cellranger/out/res_{sample_name}/outs/per_sample_outs/res_{sample_name}/vdj_b/filtered_contig_annotations.csv"))

t_contig.list <- createHTOContigList(t_contigs, 
                                     seurat_obj, 
                                     group.by = "ADT_maxID")

b_contig.list <- createHTOContigList(b_contigs, 
                                     seurat_obj, 
                                     # samples = ,
                                     group.by = "ADT_maxID")

# combineTCR - one line per cell
combined.TCR <- combineTCR(t_contig.list,
                           removeNA = FALSE, 
                           removeMulti = FALSE, 
                           filterMulti = FALSE)

head(combined.TCR[[1]])
colnames(combined.TCR[[1]])

# combineBCR - - one line per cell
# Combine using the default similarity clustering
combined.BCR <- combineBCR(b_contig.list, 
                           samples = "Patient1", 
                           filterNonproductive = TRUE, # default 
                           filterMulti = TRUE # default 
                           ) 

# The CTstrict column contains cluster IDs (e.g., "cluster.1")
head(combined.BCR.clustered[[1]][, c("barcode", "CTstrict", "IGH", "cdr3_aa1")])

# Basic Clonal Visualizations
clonalQuant(combined.TCR, 
            cloneCall="strict", 
            chain = "both", 
            scale = TRUE)


