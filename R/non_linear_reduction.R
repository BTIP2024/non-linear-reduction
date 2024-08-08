#' Non linear dimensionality reduction for scRNAseq
#' 
#' t-SNE and UMAP functions are included in this package
#' 
#' @param input is the output of RunPCA
#' @examples 
#' tsne_seurat(afterclustering.rds)
#' @export
tsne_seurat <- function(input){
   for_tsne <- readRDS(input)
   
   tsne <- Seurat::RunTSNE(for_tsne, dims = 1:10, dim.embed =2, label = TRUE)
   tsne <- Seurat::DimPlot(tsne, reduction = "tsne")
   ggplot2::ggsave(tsne, file = "tsne_seurat.png", width = 10, height = 10)
}

umap_seurat <- function(input){
   for_umap <- readRDS(input)
   
   umap <- Seurat::RunUMAP(for_umap, dims = 1:10, n.components = 2L)
   umap <- Seurat::DimPlot(umap, reduction = "umap")
   ggplot2::ggsave(umap, file = "tsne_seurat.png", width = 10, height = 10)
}