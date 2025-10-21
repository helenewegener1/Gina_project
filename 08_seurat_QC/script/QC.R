# setwd("~/ciir/people/helweg/projects/Gina_project/")

# Load libraries 
library(SeuratObject)
library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(stringr)
library(glue)
library(ggplot2)
library(DropletUtils)
library(DoubletFinder)
library(glmGamPoi)
# devtools::install_github("constantAmateur/SoupX", ref='devel')
library(SoupX)

# Load data
seurat_obj_list <- readRDS("06_seurat_load/out/seurat_obj_list.rds") # raw = prefiltering 
seurat_obj_roughQC_list <- readRDS("07_seurat_roughQC/out/seurat_obj_roughQC_list.rds")

############################# Sanity check for ADT #############################

for (sample_name in names(seurat_obj_roughQC_list)){
  
  seurat_obj_raw <- seurat_obj_roughQC_list[[sample_name]]
  
  if ("ADT" %in% names(seurat_obj_raw)){
    
    Idents(seurat_obj_raw)
    Idents(seurat_obj_raw) <- "ADT_maxID"
    RidgePlot(seurat_obj_raw, assay = "ADT", features = rownames(seurat_obj_raw[["ADT"]]))
    ggsave(glue("08_seurat_QC/plot/RidgePlot_{sample_name}.pdf"), width = 16, height = 20)
    
  }
  
}

################################################################################

# Investigate need for removal of empty droplets 
# https://bioconductor.org/packages/release/bioc/vignettes/DropletUtils/inst/doc/DropletUtils.html
for (sample_name in names(seurat_obj_list)){
  
  ######################## CHECK EMPTY DROPLETS IN PRE-ROUGH-WQ #########################
  seurat_obj_raw <- seurat_obj_list[[sample_name]]
  
  # Get count matrix
  count_mat <- seurat_obj_raw@assays$RNA@layers$counts
  
  br.out <- barcodeRanks(count_mat)
  
  # Making Barcode Rank Plot.
  pdf(glue("08_seurat_QC/plot/empty_droplets_{sample_name}_raw.pdf"), width = 8, height = 6)
  
  plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total", main=glue("{sample_name} raw (prefilter) Barcode Rank Plot"))
  o <- order(br.out$rank)
  lines(br.out$rank[o], br.out$fitted[o], col="red")
  
  abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
  abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
  legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
         legend=c("knee", "inflection"))
  
  dev.off()
  
  # Testing for empty droplets
  set.seed(100)
  
  # raw_counts <- Read10X(data.dir = glue("05_run_cellranger/out/res_{sample_name}/outs/multi/count/raw_feature_bc_matrix"))
  
  
  sc_exp_raw <- seurat_obj_raw %>% as.SingleCellExperiment()
  counts(sc_exp_raw)
  
  e.out <- emptyDrops(counts(sc_exp_raw), lower=100)
  
  # See ?emptyDrops for an explanation of why there are NA values.
  summary(e.out$FDR <= 0.001)
  
  ################### CHECK EMPTY DROPLETS IN ROUGH FILTERED ###################
  seurat_obj_roughQC <- seurat_obj_roughQC_list[[sample_name]]
  
  # Get count matrix
  count_mat <- seurat_obj_roughQC@assays$RNA@layers$counts
  
  br.out <- barcodeRanks(count_mat)
  
  # Making Barcode Rank Plot.
  pdf(glue("08_seurat_QC/plot/empty_droplets_{sample_name}_roughQC.pdf"), width = 8, height = 6)
  
  plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total", main=glue("{sample_name} roughQC-filtered Barcode Rank Plot"))
  o <- order(br.out$rank)
  lines(br.out$rank[o], br.out$fitted[o], col="red")
  
  abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
  abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
  legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
         legend=c("knee", "inflection"))
  
  dev.off()
  
}

rm(seurat_obj_list)

################################ DoubletFinder ################################ 

# Initialize final QC list 
seurat_obj_finalQC_list <- list()

for (sample_name in names(seurat_obj_roughQC_list)){
  
  # Define sample
  seurat_obj_roughQC <- seurat_obj_roughQC_list[[sample_name]]
  
  #Normalize and scale using SCTransform. Note to us: decide on vars.to.regress. 
  seurat_obj_roughQC <- SCTransform(seurat_obj_roughQC, assay = "RNA", layer = "counts", verbose = FALSE)

  #Doublet detection
  seurat_obj_roughQC <- RunPCA(seurat_obj_roughQC)
  ElbowPlot(seurat_obj_roughQC) #to determine dimentions used for following steps in doublet detection. Adjust dims. 
  seurat_obj_roughQC <- FindNeighbors(seurat_obj_roughQC, dims = 1:20)
  seurat_obj_roughQC <- FindClusters(seurat_obj_roughQC, resolution = 0.5)
  seurat_obj_roughQC <- RunUMAP(seurat_obj_roughQC, dims = 1:20)
  DimPlot(seurat_obj_roughQC)
  
  #define expected number of doublets (10x Genomics) based on the number of cells
  # Get the number of cells in the Seurat object
  num_cells <- ncol(seurat_obj_roughQC)
  print(num_cells)
  
  #identify the Pk parameter for each sample, using no ground truth strategy
  ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
  sweep.res.list <- paramSweep(seurat_obj_roughQC, PCs = 1:15, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn<- find.pK(sweep.stats)
  #visualize plot to find Pk parameter (highest peak)
  ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
    geom_point()+
    geom_line()
  
  #store max Pk as Pk variable 
  pK <- bcmvn %>% 
    filter(BCmetric == max(BCmetric)) %>%
    select(pK)
  pK <- as.numeric(as.character(pK[[1]]))
  
  #homotypic doublet proportions
  annotations <- seurat_obj_roughQC@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.061*nrow(seurat_obj_roughQC@meta.data))
  nExp.poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  #Run Doublet finder
  seurat_obj_roughQC <- doubletFinder(seurat_obj_roughQC,
                                      PCs = 1:20,
                                      pN = 0.25,
                                      pK = pK,
                                      nExp = nExp.poi.adj,
                                      # reuse.pANN = FALSE, 
                                      sct = TRUE)
  
  
  #view metadata to see single vs. double in DF.classicication
  # View(seurat_obj_roughQC@meta.data)
  #get the name of the coloumn, copy DF.classification name  
  # names(seurat_obj_roughQC@meta.data)
  
  #visualize
  DF.classification <- colnames(seurat_obj_roughQC@meta.data)[colnames(seurat_obj_roughQC@meta.data) %>% str_starts('DF')]
  
  # DimPlot(seurat_obj_roughQC, reduction = 'umap', group.by = "DF.classifications_0.25_0.3_413")
  DimPlot(seurat_obj_roughQC, reduction = 'umap', group.by = DF.classification)
  
  #number of singlet and doublet
  table(seurat_obj_roughQC@meta.data[DF.classification])
  
  # We do not do this here, but maybe later if it makes sense: Subset the Seurat object to include only "Singlet" cells
  # seurat_obj_finalQC <- seurat_obj_roughQC[, seurat_obj_roughQC@meta.data[[DF.classification]] == "Singlet"]
  # seurat_obj_finalQC_list[[sample_name]] <- seurat_obj_finalQC
  
  seurat_obj_finalQC_list[[sample_name]] <- seurat_obj_roughQC
  
}

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# N cells before and after DoubletFinder

for (sample_name in names(seurat_obj_roughQC_list)){
  print(sample_name)
  print(seurat_obj_finalQC_list[[sample_name]][[]]$DF.classification %>% table())
  print(" ")
}

############################ Ambiant RNA with SoupX ############################
# library(multtest)
# # I get the same rho for all cells... 
# sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
# 
# raw_counts <- Read10X(data.dir = glue("05_run_cellranger/out/res_{sample_name}/outs/multi/count/raw_feature_bc_matrix"))
# cell_counts <- Read10X(data.dir = glue("05_run_cellranger/out/res_{sample_name}/outs/per_sample_outs/res_{sample_name}/count/sample_filtered_feature_bc_matrix"))
# 
# sc <- SoupChannel(raw_counts$`Gene Expression`, cell_counts$`Gene Expression`, calcSoupProfile = FALSE)
# sc <-  estimateSoup(sc)
# seurat_temp <- CreateSeuratObject(counts = sc$toc)
# 
# # Perform quick clustering steps
# seurat_temp <- NormalizeData(seurat_temp, verbose = FALSE)
# seurat_temp <- FindVariableFeatures(seurat_temp, verbose = FALSE)
# seurat_temp <- ScaleData(seurat_temp, verbose = FALSE)
# seurat_temp <- RunPCA(seurat_temp, npcs = 20, verbose = FALSE)
# seurat_temp <- FindNeighbors(seurat_temp, dims = 1:20, verbose = FALSE)
# seurat_temp <- FindClusters(seurat_temp, resolution = 0.7, verbose = FALSE)
# seurat_temp <- RunUMAP(seurat_temp, dims = 1:20)
# 
# FeaturePlot(seurat_temp, features = "MKI67") 
# DimPlot(seurat_temp, group.by = "seurat_clusters")
# 
# # Extract and set the clusters
# clusters <- seurat_temp@meta.data$seurat_clusters
# sc <- setClusters(sc, clusters)
# sc <- autoEstCont(sc)
# out <- adjustCounts(sc)
# 
# # Use the 'out' matrix to create a Seurat object
# final_seurat_obj <- CreateSeuratObject(counts = out, project = "SoupX_Corrected_scRNA")
# # Add the per-cell rho to the final Seurat object
# final_seurat_obj$soupX_rho_per_cell <- sc$metaData$rho
# 
# # Visualize the contamination fraction across your UMAP
# final_seurat_obj <- NormalizeData(final_seurat_obj, verbose = FALSE)
# final_seurat_obj <- FindVariableFeatures(final_seurat_obj, verbose = FALSE)
# final_seurat_obj <- ScaleData(final_seurat_obj, verbose = FALSE)
# final_seurat_obj <- RunPCA(final_seurat_obj, npcs = 20, verbose = FALSE)
# final_seurat_obj <- FindNeighbors(final_seurat_obj, dims = 1:20, verbose = FALSE)
# final_seurat_obj <- FindClusters(final_seurat_obj, resolution = 0.5, verbose = FALSE)
# final_seurat_obj <- RunUMAP(final_seurat_obj, dims = 1:20)
# 
# FeaturePlot(final_seurat_obj, reduction = "umap", features = "soupX_rho_per_cell") 
# 
# FeaturePlot(final_seurat_obj, features = "MKI67") 




#################### Export list of filtered Seurat objects #################### 

saveRDS(seurat_obj_finalQC_list, "08_seurat_QC/out/seurat_obj_finalQC_list.rds")







