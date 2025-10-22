setwd("~/ciir/people/helweg/projects/Gina_project/")

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
library(multtest)

# Load data
seurat_obj_list <- readRDS("06_seurat_load/out/seurat_obj_list.rds") # cellranger filtered

############################# Sanity check for ADT #############################

for (sample_name in names(seurat_obj_list)){
  
  seurat_obj_raw <- seurat_obj_list[[sample_name]]
  
  if ("ADT" %in% names(seurat_obj_raw)){
    
    Idents(seurat_obj_raw)
    Idents(seurat_obj_raw) <- "ADT_maxID"
    RidgePlot(seurat_obj_raw, assay = "ADT", features = rownames(seurat_obj_raw[["ADT"]]))
    ggsave(glue("07_seurat_QC/plot/ADT_RidgePlot/RidgePlot_{sample_name}.pdf"), width = 16, height = 20)
    
  }
  
}

################################################################################

# Investigate need for removal of empty droplets 
# https://bioconductor.org/packages/release/bioc/vignettes/DropletUtils/inst/doc/DropletUtils.html
print("---------------------------------------------------------------")
print("Investigate need for removal of empty droplets")
for (sample_name in names(seurat_obj_list)){
  
  print("----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----")
  print(sample_name)
  
  ######################## CHECK EMPTY DROPLETS IN RAW #########################
  
  # Load raw (totally non-filtered) counts from cellranger 
  raw_counts <- Read10X(data.dir = glue("05_run_cellranger/out/res_{sample_name}/outs/multi/count/raw_feature_bc_matrix"))
  seurat_obj_raw <- CreateSeuratObject(counts = raw_counts)
  
  # Transform to SingleCellExperiment object
  if (length(Layers(seurat_obj_raw))){
    sc_exp_raw <- seurat_obj_raw %>% JoinLayers() %>% as.SingleCellExperiment()
  } else {
    sc_exp_raw <- as.SingleCellExperiment(seurat_obj_raw)
  }

  # Run emptyDrops
  set.seed(100)
  e.out <- emptyDrops(counts(sc_exp_raw))
  
  # See ?emptyDrops for an explanation of why there are NA values.
  is.cell <- table(e.out$FDR <= 0.01, useNA = "ifany")
  
  # Computing barcode ranks
  br.out <- barcodeRanks(counts(sc_exp_raw))
  
  # Making Barcode Rank Plot.
  pdf(glue("07_seurat_QC/plot/emptyDrops/emptyDrops_{sample_name}_raw.pdf"), width = 8, height = 6)
  
  plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total", 
       main=glue("{sample_name} raw (prefilter) Barcode Rank Plot"), 
       sub = glue("Cells: {is.cell[['TRUE']]}, Not cells: {is.cell[['FALSE']]}, NAs: {is.cell[[3]]}"))
  o <- order(br.out$rank)
  lines(br.out$rank[o], br.out$fitted[o], col="red")
  
  abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
  abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
  legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
         legend=c("knee", "inflection"))
  
  dev.off()
  
  ############################ Ambiant RNA with SoupX ############################
  
  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"

  raw_counts <- Read10X(data.dir = glue("05_run_cellranger/out/res_{sample_name}/outs/multi/count/raw_feature_bc_matrix"))
  cell_counts <- Read10X(data.dir = glue("05_run_cellranger/out/res_{sample_name}/outs/per_sample_outs/res_{sample_name}/count/sample_filtered_feature_bc_matrix"))
  
  if (is.null(names(raw_counts))){
    sc <- SoupChannel(raw_counts, cell_counts, calcSoupProfile=FALSE)
  } else if (length(names(raw_counts)) > 1) {
    sc <- SoupChannel(raw_counts$`Gene Expression`, cell_counts$`Gene Expression`, calcSoupProfile=FALSE)
  }
  
  sc <- estimateSoup(sc)

  # Perform quick clustering steps
  seurat_temp <- CreateSeuratObject(counts = sc$toc)
  seurat_temp <- NormalizeData(seurat_temp, verbose = FALSE)
  seurat_temp <- FindVariableFeatures(seurat_temp, verbose = FALSE)
  seurat_temp <- ScaleData(seurat_temp, verbose = FALSE)
  seurat_temp <- RunPCA(seurat_temp, npcs = 20, verbose = FALSE)
  seurat_temp <- FindNeighbors(seurat_temp, dims = 1:20, verbose = FALSE)
  seurat_temp <- FindClusters(seurat_temp, resolution = 0.7, verbose = FALSE)
  seurat_temp <- RunUMAP(seurat_temp, dims = 1:20)
  
  DimPlot(seurat_temp, group.by = "seurat_clusters") + labs(title = "Pre-SoupX")
  ggsave(glue("07_seurat_QC/plot/SoupX/pre_SoupX_{sample_name}.pdf"), width = 7, height = 6)

  # Extract and set the clusters
  clusters <- seurat_temp@meta.data$seurat_clusters
  sc <- setClusters(sc, clusters)
  sc <- autoEstCont(sc)
  out <- adjustCounts(sc)
  
  # Get the UMAP coordinates and add them to the SoupChannel object
  # umap_coords <- seurat_temp@reductions$umap@cell.embeddings
  # sc <- setDR(sc, umap_coords)
  # plotMarkerMap(sc, "MKI67")

  # Use the 'out' matrix to create a Seurat object
  final_seurat_obj <- CreateSeuratObject(counts = out, project = "SoupX_Corrected_scRNA")
  # Add the per-cell rho to the final Seurat object
  final_seurat_obj$soupX_rho_per_cell <- sc$metaData$rho

  # Visualize the contamination fraction across your UMAP
  final_seurat_obj <- NormalizeData(final_seurat_obj, verbose = FALSE)
  final_seurat_obj <- FindVariableFeatures(final_seurat_obj, verbose = FALSE)
  final_seurat_obj <- ScaleData(final_seurat_obj, verbose = FALSE)
  final_seurat_obj <- RunPCA(final_seurat_obj, npcs = 20, verbose = FALSE)
  final_seurat_obj <- FindNeighbors(final_seurat_obj, dims = 1:20, verbose = FALSE)
  final_seurat_obj <- FindClusters(final_seurat_obj, resolution = 0.7, verbose = FALSE)
  final_seurat_obj <- RunUMAP(final_seurat_obj, dims = 1:20)

  DimPlot(final_seurat_obj, group.by = "seurat_clusters") + labs(title = "Post-SoupX")
  ggsave(glue("07_seurat_QC/plot/SoupX/post_SoupX_{sample_name}.pdf"), width = 7, height = 6)

}

################################ DoubletFinder on cellranger filtered ################################ 

# Initialize final QC list 
seurat_obj_DoubletFinder <- list()

print("---------------------------------------------------------------")
print("DoubletFinder")

for (sample_name in names(seurat_obj_list)){
  
  print("----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----")
  print(sample_name)
  
  # Define sample
  seurat_obj <- seurat_obj_list[[sample_name]]
  
  # Seurat workflow
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  #Doublet detection
  seurat_obj <- RunPCA(seurat_obj)
  ElbowPlot(seurat_obj) #to determine dimentions used for following steps in doublet detection. Adjust dims. 
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
  DimPlot(seurat_obj)
  
  #define expected number of doublets (10x Genomics) based on the number of cells
  # Get the number of cells in the Seurat object
  num_cells <- ncol(seurat_obj)
  print(num_cells)
  
  #identify the Pk parameter for each sample, using no ground truth strategy
  ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
  sweep.res.list <- paramSweep(seurat_obj, PCs = 1:20, sct = TRUE)
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
  annotations <- seurat_obj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.061*nrow(seurat_obj@meta.data))
  nExp.poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  #Run Doublet finder
  seurat_obj <- doubletFinder(seurat_obj,
                              PCs = 1:20,
                              pN = 0.25,
                              pK = pK,
                              nExp = nExp.poi.adj,
                              # reuse.pANN = FALSE, 
                              sct = TRUE)
  
  
  #view metadata to see single vs. double in DF.classicication
  # View(seurat_obj@meta.data)
  #get the name of the coloumn, copy DF.classification name  
  # names(seurat_obj@meta.data)
  
  # visualize
  DF.classification <- colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) %>% str_starts('DF')]
  
  # Number of singlet and doublet - Add to plot
  result <- table(seurat_obj@meta.data[DF.classification], useNA = "ifany")
  
  # Plot
  DimPlot(seurat_obj, reduction = 'umap', group.by = DF.classification) + 
    labs(title = "DoubletFinder", subtitle = glue("N doublets: {result[[1]]}, N singlets: {result[[2]]}"))
  ggsave(glue("07_seurat_QC/plot/DoubletFinder/DoubletFinder_{sample_name}.pdf"), width = 7, height = 6)
  
  FeaturePlot(seurat_obj, reduction = 'umap', features = "MKI67") + 
    labs(title = "MKI67", subtitle = "MKI67 is a proliferation marker")
  ggsave(glue("07_seurat_QC/plot/DoubletFinder/MKI67_{sample_name}.pdf"), width = 7, height = 6)
  
  
  # We do not do this here, but maybe later if it makes sense: Subset the Seurat object to include only "Singlet" cells
  # seurat_obj_finalQC <- seurat_obj[, seurat_obj@meta.data[[DF.classification]] == "Singlet"]
  # seurat_obj_finalQC_list[[sample_name]] <- seurat_obj_finalQC
  
  seurat_obj_DoubletFinder[[sample_name]] <- seurat_obj

}

#################### Export list of Seurat objects with QC metrices in metadata #################### 

saveRDS(seurat_obj_DoubletFinder, "07_seurat_QC/out/seurat_obj_QC_metrics.rds")







