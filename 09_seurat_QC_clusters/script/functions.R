# Function to do scatter plot of markers for one lineage (e.g. T cells) VS another (e.g. B cells)
# colored by cell cycle phase and doublet classification (singlet/doublet)
# Idea is to detect cells that are marked as doublets but are in fact proliferating cells that we want to save. 
# Also, we are more sure that it is a doublet if it expresses markers from two different lineages 
# (the droplet might contain a T cell and a B cell)
# It is difficult to detect a droplet containing two T cells... 
doublet_dual_lineages <- function(seurat_obj, sample_name, marker_1, marker_2) {
  
  seurat_obj[[]] %>% 
    ggplot(aes(x = !!sym(marker_1), 
               y = !!sym(marker_2), 
               color = scDblFinder.class,
               shape = Phase)) +
    geom_point(alpha = 0.5) + 
    scale_color_manual(values = c("grey", "red")) + 
    theme_bw()
  
  ggsave(glue("09_seurat_QC_clusters/plot/{sample_name}/doublet/{sample_name}_doublet_Phase_{marker_1}_VS_{marker_2}_red.png"))
  
  seurat_obj[[]] %>% 
    mutate(doublet_Phase = glue("{scDblFinder.class}_{Phase}")) %>% 
    ggplot(aes(x = !!sym(marker_1), 
               y = !!sym(marker_2), 
               color = doublet_Phase,
               shape = scDblFinder.class)) +
    geom_point(alpha = 0.5) + 
    scale_color_manual(values = c("grey", "blue", "green", "grey","grey","grey")) +
    theme_bw()
  
  ggsave(glue("09_seurat_QC_clusters/plot/{sample_name}/doublet/{sample_name}_doublet_Phase_{marker_1}_VS_{marker_2}.png"))
  
}

# Function to do violin plots of nFeature_RNA and nCount_RNA
# split on x axis by doublet/singlet class
# color by cell cycle phase
# Idea is to detect cells that are marked as doublets but are infact proliferating cells that we want to save. 
doublet_N_genes <- function(seurat_obj, sample_name) {
  
  # nFeature_RNA
  seurat_obj[[]] %>% 
    ggplot(aes(x = scDblFinder.class, 
               y = nFeature_RNA)) + 
    geom_violin() + 
    geom_jitter(aes(color = Phase), width = 0.3, alpha = 0.3) + 
    scale_color_manual(values = c("grey", "red", "blue")) + 
    theme_bw()
  
  ggsave(glue("09_seurat_QC_clusters/plot/{sample_name}/doublet/{sample_name}_doublet_VS_nFeature_RNA.png"))
  
  seurat_obj[[]] %>% 
    ggplot(aes(x = scDblFinder.class, 
               y = nCount_RNA)) + 
    geom_violin() + 
    geom_jitter(aes(color = Phase), width = 0.3, alpha = 0.3) + 
    scale_color_manual(values = c("grey", "red", "blue")) + 
    theme_bw()
  
  ggsave(glue("09_seurat_QC_clusters/plot/{sample_name}/doublet/{sample_name}_doublet_VS_nCount_RNA.png"))
  
}
