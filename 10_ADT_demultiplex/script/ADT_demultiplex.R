getwd()

library(SeuratObject)
library(DropletUtils)
library(Seurat)
library(tidyverse)
library(glue)
library(patchwork)
library(readxl)
# library(ztable)
library(pheatmap)

# Load data
seurat_obj_nonDC_list <- readRDS("09_annotation_pre_integration/out/seurat_obj_nonDC_list.rds")
seurat_obj <- seurat_obj_nonDC_list[["HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"]]
# Gina_seurat_obj <- readRDS("00_data/Gina_HH117_PP_broadAnn.rds")

# CLR normalization of ADT
seurat_obj <- NormalizeData(seurat_obj, assay = "ADT", normalization.method = "CLR")

# Demultiplexing with HTODemux and MULTIseqDemux
seurat_obj <- HTODemux(seurat_obj, assay = "ADT", positive.quantile = 0.99)
seurat_obj <- MULTIseqDemux(seurat_obj, assay = "ADT", autoThresh = TRUE)

# seurat_obj$MULTI_ID
# seurat_obj$ADT_classification 

# Mutate columns so result from HTODemux and MULTIseqDemux become more comparable 
seurat_obj@meta.data$ADT_ID <- lapply(seurat_obj@meta.data$ADT_classification, function(x) {ifelse(str_detect(x, "_"), "Doublet", x)}) %>% unlist()
seurat_obj$MULTI_classification.global <- seurat_obj$MULTI_ID %>% str_replace("Fol-\\d+", "Singlet")

# Look at results
seurat_obj$ADT_ID %>% table()
seurat_obj$MULTI_ID %>% table()

seurat_obj$ADT_classification.global %>% table()
seurat_obj$MULTI_classification.global %>% table()

##################################### Plot ##################################### 

# Plot Heatmap
HTOHeatmap(seurat_obj, assay = "ADT", ncells = 5000)

# Plot density of CLR values for each ADT
# Get CLR-normalized data
adt_data <- GetAssayData(seurat_obj, assay = "ADT", layer = "data")
dim(adt_data)

# Compute mean CLR value per hashtag
adt_means <- Matrix::rowMeans(adt_data)

# Summarize in a table
adt_summary <- data.frame(
  ADT = names(adt_means),
  mean_CLR = adt_means,
  median_CLR = apply(adt_data, 1, median),
  max_CLR = apply(adt_data, 1, max)
) %>%
  arrange(desc(mean_CLR))

adt_summary %>% arrange(ADT)

# Weak tags are the ones with low CLR mean 
# weak_tags <- adt_summary %>% arrange(mean_CLR) %>% head(6) %>% pull(ADT)
# RidgePlot(Gina_seurat_obj_other, assay = "ADT", features = weak_tags, ncol = 3)

# log count value range
seurat_obj[["ADT"]]$counts %>% t() %>% as.data.frame() %>% 
  pivot_longer(cols = starts_with("Fol"), names_to = "ADT", values_to = "counts") %>% 
  mutate(ADT_CLR_mean = glue('{ADT} CLR_mean: {round(adt_summary[ADT,"mean_CLR"], 2)}')) %>% 
  ggplot(aes(x = log(counts), fill = ADT_CLR_mean)) + 
  geom_density(alpha = 0.5) + 
  facet_wrap(vars(ADT_CLR_mean)) + 
  theme_bw() + 
  theme(legend.position = "none") 

ggsave(glue("10_ADT_demultiplex/plot/ADT_log_counts_density.png"), width = 10, height = 7)

# CLR value range
seurat_obj[["ADT"]]$data %>% t() %>% as.data.frame() %>% 
  pivot_longer(cols = starts_with("Fol"), names_to = "ADT", values_to = "CLR") %>% 
  mutate(ADT_CLR_mean = glue('{ADT} CLR_mean: {round(adt_summary[ADT,"mean_CLR"], 2)}')) %>% 
  ggplot(aes(x = CLR, fill = ADT_CLR_mean)) + 
  geom_density(alpha = 0.5) + 
  facet_wrap(vars(ADT_CLR_mean)) + 
  theme_bw() + 
  theme(legend.position = "none") 

ggsave(glue("10_ADT_demultiplex/plot/ADT_CLR_values_density.png"), width = 10, height = 7)

# log CLR range
seurat_obj[["ADT"]]$data %>% t() %>% as.data.frame() %>% 
  pivot_longer(cols = starts_with("Fol"), names_to = "ADT", values_to = "CLR") %>% 
  mutate(ADT_CLR_mean = glue('{ADT} CLR_mean: {round(adt_summary[ADT,"mean_CLR"], 2)}')) %>% 
  ggplot(aes(x = log(CLR), fill = ADT_CLR_mean)) + 
  geom_density(alpha = 0.5) + 
  facet_wrap(vars(ADT_CLR_mean)) + 
  theme_bw() + 
  theme(legend.position = "none") 

ggsave(glue("10_ADT_demultiplex/plot/ADT_log_CLR_values_density.png"), width = 10, height = 7)

# 99th quantile (default thershold for HTODemux)
adt_long <- seurat_obj[["ADT"]]$data %>% t() %>% as.data.frame() %>% 
  rownames_to_column("cell") %>% 
  pivot_longer(cols = starts_with("Fol"), names_to = "ADT", values_to = "CLR") %>% 
  mutate(ADT_CLR_mean = glue('{ADT} CLR_mean: {round(adt_summary[ADT,"mean_CLR"], 2)}'))

# Compute 99th quantile per ADT
quantiles <- adt_long %>%
  group_by(ADT) %>%
  summarize(q99 = quantile(CLR, 0.99))

# Bool of cells above the 99th quantile 
n_above_quantile <- adt_long %>% left_join(quantiles, by = "ADT") %>% ungroup()
n_above_quantile <- n_above_quantile %>% mutate(is_q99 = ifelse(CLR >= q99, TRUE, FALSE))

# N cells above the 99th quantile 
n_q99 <- (ncol(seurat_obj) / 100) * (100 - 99) 
n_q99 <- ceiling(n_q99)

# ADT range with vertical line at 99th quantile 
n_above_quantile %>% 
  ggplot(aes(x = log(CLR), fill = ADT_CLR_mean)) + 
  geom_density(alpha = 0.5) + 
  geom_vline(aes(xintercept = log(q99)), linetype = "dashed") +
  facet_wrap(vars(ADT_CLR_mean)) +
  theme_bw() + 
  theme(legend.position = "none") + 
  labs(subtitle = glue("q99, n = {n_q99}"))

ggsave(glue("10_ADT_demultiplex/plot/ADT_log_CLR_values_q99_density.png"), width = 10, height = 7)

# # FeatureScatter Plots 
# for (i in 1:(n-1)) {
#   for (j in (i+1):n) {
#     fol1 <- fols[i]
#     fol2 <- fols[j]
#     
#     p <- FeatureScatter(seurat_obj, feature1 = fol1, feature2 = fol2, slot = "data") # data == CLR transformed
#     
#     # Save plot
#     ggsave(
#       filename = glue("10_ADT_demultiplex/plot/FeatureScatter/FeatureScatter_{fol1}_{fol2}.png"),
#       plot = p,
#       width = 10,
#       height = 8
#     )
#   }
# }

###################### CONSENSUS OF HTODemux AND HTODemux ######################

seurat_obj$ADT_ID

# Compare result of methods 
table(seurat_obj$ADT_ID, seurat_obj$MULTI_ID)

table(seurat_obj$ADT_ID, seurat_obj$MULTI_ID) %>% 
  heatmap(Colv = NA, Rowv = NA, xlab = "HTODemux", ylab = "MULTIseqDemux")
ggsave(glue("10_ADT_demultiplex/plot/compare_methods_labels_heatmap.png"), width = 8, height = 7)

# Consensus data used to plot and extract consensus column for final seurat object 
df_consens <- seurat_obj@meta.data %>% 
  select(ADT_classification, ADT_classification.global, ADT_ID, ADT_maxID, ADT_secondID, MULTI_ID, MULTI_classification, MULTI_classification.global) %>% 
  # Extract MULTI_maxID and MULTI_secondID to compare to HTODemux annotations 
  mutate(
    MULTI_maxID = str_split_i(MULTI_classification, "_", 1),
    MULTI_secondID = str_split_i(MULTI_classification, "_", 2)
  ) %>% 
  mutate(
    
    # Everything same, also same exact doublets, else NA
    ADT_consensus_very_hard = ifelse(ADT_classification == MULTI_classification, ADT_classification, NA),
    
    # Same singlet, negative and doublet detection, else NA
    ADT_consensus_hard = ifelse(ADT_ID == MULTI_ID, ADT_ID, NA),
    
    # Since HTODemux detects many doublets and MULTIseqDemux detects more negatives. 
    ADT_consensus_medium = case_when(
      # if HTODemux and MULTIseqDemux agree - keep label
      ADT_ID == MULTI_ID ~ ADT_ID,
      
      # if HTODemux detects a singlet and MULTIseqDemux detects a singlet - let it be a negative if the top follicles are NOT identical
      ADT_classification.global == "Singlet" & MULTI_classification.global == "Singlet" & ADT_classification != MULTI_maxID ~ "Negative",
      
      # if HTODemux detects a doublet and MULTIseqDemux detects a singlet - let it be a singlet IF ADT_maxID == MULTI_ID
      ADT_classification.global == "Doublet" & MULTI_classification.global == "Singlet" & ADT_maxID == MULTI_ID ~ MULTI_ID,
      # else if ADT_maxID != MULTI_ID
      ADT_classification.global == "Doublet" & MULTI_classification.global == "Singlet" & ADT_maxID != MULTI_ID ~ "Doublet",
      
      # if HTODemux detects a singlet and MULTIseqDemux detects a negative - let it be a singlet IF ADT_classification == MULTI_maxID
      # MULTI_maxID will be Negative if MULTI_ID is Negative 
      # ADT_classification.global == "Singlet" & MULTI_ID == "Negative" & ADT_classification == MULTI_maxID ~ ADT_classification,
      # else if ADT_maxID != MULTI_ID
      # ADT_classification.global == "Singlet" & MULTI_ID == "Negative" ~ "Negative",
      ADT_classification.global == "Singlet" & MULTI_ID == "Negative" ~ ADT_classification,
      
      # if HTODemux detects a doublet and MULTIseqDemux detects a negative - let it be a negative
      ADT_classification.global == "Doublet" & MULTI_ID == "Negative" ~ "Negative",
      
      # if HTODemux detects a doublet and MULTIseqDemux detects a doublet - let it be a singlet if the two first follicles are identical but the second ones aren't.
      ADT_classification.global == "Doublet" & MULTI_ID == "Doublet" & ADT_maxID == MULTI_maxID & ADT_secondID != MULTI_secondID ~ ADT_maxID,
      # if HTODemux detects a doublet and MULTIseqDemux detects a doublet - let it be a doublet if the two first and the second follicles are identical.
      ADT_classification.global == "Doublet" & MULTI_ID == "Doublet" & ADT_maxID == MULTI_maxID & ADT_secondID == MULTI_secondID ~ "Doublet",
      # if HTODemux detects a doublet and MULTIseqDemux detects a doublet - let it be a doublet
      ADT_classification.global == "Doublet" & MULTI_ID == "Doublet" & ADT_maxID != MULTI_maxID ~ "Doublet"

    ),
    ADT_consensus_soft = case_when(
      # if HTODemux and MULTIseqDemux agree - keep label
      ADT_ID == MULTI_ID ~ ADT_ID,
      
      # if HTODemux detects a singlet and MULTIseqDemux detects a singlet - let it be a negative if the top follicles are NOT identical
      ADT_classification.global == "Singlet" & MULTI_classification.global == "Singlet" & ADT_classification != MULTI_maxID ~ "Negative",
      
      # if HTODemux detects a doublet and MULTIseqDemux detects a singlet - let it be a singlet
      ADT_classification.global == "Doublet" & MULTI_classification.global == "Singlet"  ~ MULTI_ID,
      # if HTODemux detects a singlet and MULTIseqDemux detects a negative - let it be a singlet
      ADT_classification.global == "Singlet" & MULTI_ID == "Negative" ~ ADT_classification,
      
      # if HTODemux detects a doublet and MULTIseqDemux detects a negative - let it be a negative
      ADT_classification.global == "Doublet" & MULTI_ID == "Negative" ~ "Negative",
      
      # if HTODemux detects a doublet and MULTIseqDemux detects a doublet - let it be a singlet if the two first follicles are identical.
      ADT_classification.global == "Doublet" & MULTI_ID == "Doublet" & ADT_maxID == MULTI_maxID  ~ ADT_maxID,
      # if HTODemux detects a doublet and MULTIseqDemux detects a doublet - let it be a doublet
      ADT_classification.global == "Doublet" & MULTI_ID == "Doublet" & ADT_maxID != MULTI_maxID ~ "Doublet"
      
    )
  )

head(df_consens) 

# Tables
table(df_consens$ADT_consensus_very_hard, useNA = "always") # same doublet 
table(df_consens$ADT_consensus_hard, useNA = "always")
table(df_consens$ADT_consensus_medium, useNA = "always")
table(df_consens$ADT_consensus_soft, useNA = "always")

df_consens$ADT_classification.global %>% table()
df_consens$MULTI_classification.global %>% table()

df_consens$ADT_consensus_hard.global <- df_consens$ADT_consensus_hard %>% str_replace("Fol-\\d+", "Singlet") 
df_consens$ADT_consensus_medium.global <- df_consens$ADT_consensus_medium %>% str_replace("Fol-\\d+", "Singlet") 
df_consens$ADT_consensus_soft.global <- df_consens$ADT_consensus_soft %>% str_replace("Fol-\\d+", "Singlet")

table(df_consens$ADT_consensus_hard.global, useNA = "always")
table(df_consens$ADT_consensus_medium.global, useNA = "always")
table(df_consens$ADT_consensus_soft.global, useNA = "always")

# When MULTIseqDemux detects Negative, what does HTODemux detect?
df_consens %>% 
  filter(MULTI_ID == "Negative") %>% 
  ggplot(aes(x = ADT_ID)) + 
  geom_bar() +
  theme_bw() + 
  labs(title = "MULTIseqDemux = Negative", 
       subtitle = "When MULTIseqDemux detects Negative, what does HTODemux detect?")

ggsave("10_ADT_demultiplex/plot/MULTIseqDemux_Negative_barplot.png", width = 10, height = 7)

# When HTODemux detects Doublet, what is the max and second tag?
df_consens %>% 
  filter(ADT_classification.global == "Doublet") %>% 
  ggplot(aes(x = ADT_maxID, 
             color = ADT_secondID)) + 
  geom_bar() +
  theme_bw() + 
  labs(title = "HTODemux = Doublet", 
       subtitle = "When HTODemux detects Doublet, what is the max and second tag?")

ggsave("10_ADT_demultiplex/plot/HTODemux_Doublet.png", width = 10, height = 7)

#################### Add consensus colums to seurat object #####################

seurat_obj@meta.data %>% dim()
df_consens %>% dim()

seurat_obj@meta.data$ADT_consensus_very_hard <- df_consens$ADT_consensus_very_hard
seurat_obj@meta.data$ADT_consensus_hard <- df_consens$ADT_consensus_hard
seurat_obj@meta.data$ADT_consensus_medium <- df_consens$ADT_consensus_medium
seurat_obj@meta.data$ADT_consensus_soft<- df_consens$ADT_consensus_soft

seurat_obj@meta.data$ADT_consensus_hard.global<- df_consens$ADT_consensus_hard.global
seurat_obj@meta.data$ADT_consensus_medium.global<- df_consens$ADT_consensus_medium.global
seurat_obj@meta.data$ADT_consensus_soft.global <- df_consens$ADT_consensus_soft.global

############################### ADT Ridge Plots ################################ 

idents <- c("ADT_maxID", "ADT_ID", "ADT_consensus_medium", "ADT_consensus_soft")
  
for (ident in idents){
  
  Idents(seurat_obj) <- ident
  RidgePlot(seurat_obj, assay = "ADT", features = rownames(seurat_obj[["ADT"]])) 
  ggsave(glue("10_ADT_demultiplex/plot/RidgePlot_{ident}.pdf"), width = 16, height = 20)
  
}

################################### DimPlot ####################################

DimPlot(seurat_obj, group.by = "ADT_classification.global")
DimPlot(seurat_obj, group.by = "MULTI_classification.global")

# Consensus
DimPlot(seurat_obj, group.by = "ADT_consensus_hard.global")
DimPlot(seurat_obj, group.by = "ADT_consensus_medium.global")
DimPlot(seurat_obj, group.by = "ADT_consensus_soft.global")

seurat_obj$ADT_consensus_soft

################################### tSNE ####################################

Idents(seurat_obj) <- "ADT_consensus_soft.global"

# First, we will remove negative cells from the object
seurat_obj.subset <- subset(seurat_obj, idents = "Negative", invert = TRUE)

# Calculate a tSNE embedding of the HTO data
DefaultAssay(seurat_obj.subset) <- "ADT"
seurat_obj.subset <- ScaleData(seurat_obj.subset, features = rownames(seurat_obj.subset), verbose = FALSE)
seurat_obj.subset <- RunPCA(seurat_obj.subset, features = rownames(seurat_obj.subset), approx = FALSE)
seurat_obj.subset <- RunUMAP(seurat_obj.subset, dims = 1:8, perplexity = 100)
seurat_obj.subset <- RunTSNE(seurat_obj.subset, dims = 1:8, perplexity = 100)

DimPlot(seurat_obj.subset, reduction = "umap")
ggsave("10_ADT_demultiplex/plot/DimPlot_Doublet_Singlet_umap.png", width = 9, height = 7)

DimPlot(seurat_obj.subset, reduction = "tsne")
ggsave("10_ADT_demultiplex/plot/DimPlot_Doublet_Singlet_tsne.png", width = 9, height = 7)

####################### Export ADT demultiplexed objects ####################### 

saveRDS(seurat_obj, "10_ADT_demultiplex/out/seurat_obj_ADT_demultiplexed.rds")





