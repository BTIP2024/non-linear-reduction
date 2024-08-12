#' Non linear dimensionality reduction for scRNAseq
#' 
#' t-SNE and UMAP functions are included in this package
#' 
#' @param input is the output of RunPCA
#' @examples 
#' tsne_seurat(afterclustering.rds)
#' @export
clusters_seurat <- function(input){
   clustering <- readRDS(input)
   clustering <- Seurat::RunPCA(clustering, features = Seurat::VariableFeatures(object = clustering))
   image4 <- Seurat::FindNeighbors(clustering, dims= 1:15)
   image4 <- Seurat::FindClusters(image4, resolution = c(0.1, 0.3, 0.5, 0.7, 1))
   saveRDS(image4, file = "clustersperresolution.rds")
   
   image4 <- Seurat::DimPlot(image4, group.by = "RNA_snn_res.0.3", label = TRUE)
   
   ggplot2::ggsave(image4, file = "dimplot.png", width = 12, height = 10)
}

tsne_seurat <- function(input){
   for_tsne <- readRDS(input)
   for_3d <- for_tsne
   
   tsne <- Seurat::RunTSNE(for_tsne, dims = 1:10, dim.embed =2, label = TRUE)
   tsne <- Seurat::DimPlot(tsne, reduction = "tsne")
   ggplot2::ggsave(tsne, file = "tsne_seurat.png", width = 10, height = 10)
   
   
   #for 3D plots
   tsne_3d <- Seurat::RunTSNE(for_3d, dims = 1:10, dim.embed = 3)
   
   tsne_1 <- tsne_3d[["tsne"]]@cell.embeddings[,1]
   tsne_2 <- tsne_3d[["tsne"]]@cell.embeddings[,2]
   tsne_3 <- tsne_3d[["tsne"]]@cell.embeddings[,3]
   
   plot.data <- Seurat::FetchData(object = tsne_3d, vars = c("tSNE_1", "tSNE_2", "tSNE_3", "seurat_clusters"))
   
   plot.data$label <- paste(rownames(plot.data))
   
   plotin3d <- plotly::plot_ly(data = plot.data, 
                               x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, 
                               color = ~seurat_clusters,
                               type = "scatter3d", 
                               mode = "markers", 
                               marker = list(size = 5, width=2), 
                               text=~label, 
                               hoverinfo="text")
   
   htmltools::save_html(plotin3d, file = "tsne_3dplot.html")
   
}

umap_seurat <- function(input){
   for_umap <- readRDS(input)
   
   umap <- Seurat::RunUMAP(for_umap, dims = 1:10, n.components = 3L)
   umap <- Seurat::DimPlot(umap, reduction = "umap")
   ggplot2::ggsave(umap, file = "tsne_seurat.png", width = 10, height = 10)
   
   plot.data <- Seurat::FetchData(object = for_umap, vars = c("umap_1", "umap_2", "umap_3", "seurat_clusters"))
   
   plot.data$label <- paste(rownames(plot.data))
   
   fig <- plot_ly(data = plot.data, 
                  x = ~umap_1, y = ~umap_2, z = ~umap_3, 
                  color = ~seurat_clusters,
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 5, width=2), # controls size of points
                  text=~label, 
                  hoverinfo="text")
   
   htmltools::save_html(fig, file = "umap_3dplot.html")
}