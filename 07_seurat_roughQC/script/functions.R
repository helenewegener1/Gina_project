
plot_qc <- function(seurat_obj, sample_name, n_cells, version = "raw", filtering = ""){
  
  p_ribo <- VlnPlot(seurat_obj, features = "percent.ribo", layer = "counts") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = 'none')
  p_mt <- VlnPlot(seurat_obj, features = "percent.mt", layer = "counts") + geom_hline(yintercept = 5, color = "black") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = 'none')
  p_feature <- VlnPlot(seurat_obj, features = "nFeature_RNA", layer = "counts") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) + theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = 'none')
  p_count <- VlnPlot(seurat_obj, features = "nCount_RNA", layer = "counts") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) + theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = 'none')
  
  p_final <- (p_ribo | p_mt) / (p_feature | p_count) + 
    plot_annotation(title = glue("{tools::toTitleCase(version)} QC plots of sample {sample_name}"),
                    caption = glue("N cells: {n_cells}"))
  
  if (filtering != ""){
    p_final <- p_final + 
      plot_annotation(subtitle = filtering)
  }
  
  ggsave(plot = p_final,
         filename = glue("07_seurat_roughQC/plot/{version}_QC_plot_{sample_name}.png"), 
         width = 9, 
         height = 8)
  
  return(p_final)
  
}
