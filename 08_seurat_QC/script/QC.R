# setwd("~/ciir/people/helweg/projects/Gina_project/")

# Load libraries 
library(SeuratObject)
library(Seurat)
library(dplyr)
library(stringr)
library(glue)
library(ggplot2)
library(DropletUtils)
library(DoubletFinder)
library(glmGamPoi)

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
  
  ######################## CHECK EMPTY DROPLETS IN RAW #########################
  seurat_obj_raw <- seurat_obj_list[[sample_name]]
  
  # Get count matrix
  count_mat <- seurat_obj_raw@assays$RNA@layers$counts
  
  br.out <- barcodeRanks(count_mat)
  
  # Making Barcode Rank Plot.
  pdf(glue("08_seurat_QC/plot/{sample_name}_empty_droplets_plot_raw.pdf"), width = 8, height = 6)
  
  plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total", main=glue("{sample_name} raw (prefilter) Barcode Rank Plot"))
  o <- order(br.out$rank)
  lines(br.out$rank[o], br.out$fitted[o], col="red")
  
  abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
  abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
  legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
         legend=c("knee", "inflection"))
  
  dev.off()
  
  ################### CHECK EMPTY DROPLETS IN ROUGH FILTERED ###################
  seurat_obj_roughQC <- seurat_obj_roughQC_list[[sample_name]]
  
  # Get count matrix
  count_mat <- seurat_obj_roughQC@assays$RNA@layers$counts
  
  br.out <- barcodeRanks(count_mat)
  
  # Making Barcode Rank Plot.
  pdf(glue("08_seurat_QC/plot/{sample_name}_empty_droplets_plot_roughQC.pdf"), width = 8, height = 6)
  
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

################################ Doublet finder ################################ 

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
  
  # Subset the Seurat object to include only "Singlet" cells
  seurat_obj_finalQC <- seurat_obj_roughQC[, seurat_obj_roughQC@meta.data[[DF.classification]] == "Singlet"]
  
  seurat_obj_finalQC_list[[sample_name]] <- seurat_obj_finalQC
  
}

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# N cells before and after DoubletFinder

for (sample_name in names(seurat_obj_roughQC_list)){
  print(sample_name)
  print(seurat_obj_roughQC_list[[sample_name]] %>% ncol()) 
  print(seurat_obj_finalQC_list[[sample_name]] %>% ncol())
  print(" ")
}

#################### Export list of filtered Seurat objects #################### 

saveRDS(seurat_obj_finalQC_list, "08_seurat_QC/out/seurat_obj_finalQC_list.rds")







