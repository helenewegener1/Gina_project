getwd()

library(SeuratObject)
library(Seurat)
library(dplyr)
library(stringr)
library(tibble)
library(tidyr)
library(glue)
library(ggplot2)
library(patchwork)
library(readxl)
library(ztable)
library(pheatmap)

# Load data
# seurat_obj_list <- readRDS("07_seurat_QC/out/seurat_obj_QC.rds")
seurat_obj_list <- readRDS("08_seurat_QC_filtering/out/seurat_obj_roughQC_list.rds")
seurat_obj <- seurat_obj_list[["HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"]]

Gina_seurat_obj <- readRDS("00_data/Gina_HH117_PP_broadAnn.rds")

# Explore data
dim(seurat_obj@assays$RNA)
dim(Gina_seurat_obj@assays$RNA)

DimPlot(Gina_seurat_obj, group.by = )

# Seurat workflow
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, npcs = 20, verbose = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.1, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

DimPlot(seurat_obj)

# Clusters
Gina_seurat_obj@meta.data$Celltype 
Gina_seurat_obj


# Gina seurat object 

Gina_seurat_obj@meta.data$Celltype
Gina_seurat_obj@meta.data$ADT_maxID
Gina_seurat_obj@meta.data$ADT_classification
Gina_seurat_obj@meta.data$ADT_classification.global

x_1 <- table(Gina_seurat_obj@meta.data$Celltype, Gina_seurat_obj@meta.data$ADT_classification.global)
x_2 <- table(Gina_seurat_obj@meta.data$Celltype, Gina_seurat_obj@meta.data$ADT_classification)
x_3 <- table(Gina_seurat_obj@meta.data$Celltype, Gina_seurat_obj@meta.data$ADT_maxID)

heatmap(x_1, Colv = NA, Rowv = NA)
# heatmap(x_2, Colv = NA, Rowv = NA)
heatmap(x_3, Colv = NA, Rowv = NA)

# Plot
x_1 %>% t() %>% as.data.frame.matrix() %>% rownames_to_column("global_class") %>% 
  pivot_longer(cols = c("Mem.B cells", "Tfh cells", "GC.B cells", "Naive.B cells", "PBs", "DCs")) %>% 
  ggplot(aes(y = name, x = value, fill = global_class)) + 
  geom_col(position = "dodge") + 
  theme_bw() + 
  scale_fill_manual(values = c("blue4", "red3", "lightblue")) + 
  labs(x = "Count", y = "")
ggsave("07_seurat_QC/ADT_explore/global_class.png", width = 12, height = 8)

# x_2 %>% t() %>% as.data.frame.matrix() %>% rownames_to_column("class") %>% 
#   # pivot_longer(cols = c("Mem.B cells", "Tfh cells", "GC.B cells", "Naive.B cells", "PBs", "DCs")) %>% 
#   ggplot(aes(x = DCs, y = class)) + 
#   geom_col(position = "dodge") + 
#   theme_bw()  
#   # scale_fill_manual(values = c("blue4", "red3", "lightblue"))

x_3 %>% t() %>% as.data.frame.matrix() %>% rownames_to_column("maxID") %>% 
  ggplot(aes(x = DCs, y = maxID)) + 
  geom_col(position = "dodge") + 
  theme_bw() + 
  labs(x = "Count")
ggsave("07_seurat_QC/ADT_explore/maxID_DCs.png", width = 12, height = 8)

x_3 %>% t() %>% as.data.frame.matrix() %>% rownames_to_column("maxID") %>% 
  pivot_longer(cols = c("Mem.B cells", "Tfh cells", "GC.B cells", "Naive.B cells", "PBs", "DCs")) %>%
  ggplot(aes(y = name, x = value, fill = maxID)) + 
  geom_col(position = "dodge", ) + 
  theme_bw() + 
  labs(x = "Count", y = "")
ggsave("07_seurat_QC/ADT_explore/maxID_all_cells.png", width = 12, height = 8)

# Doublets
df_doublets <- Gina_seurat_obj@meta.data %>% 
  filter(ADT_classification.global == "Doublet") %>% 
  select(starts_with("ADT"))

table(df_doublets$ADT_maxID, df_doublets$ADT_secondID) %>% 
  heatmap(Colv = NA, Rowv = NA, main = "ADT_maxID vs ADT_secondID for Doublets", xlab = "ADT_maxID", ylab = "ADT_secondID")
ggsave("07_seurat_QC/ADT_explore/doublets_heatmap.png", width = 12, height = 8)






