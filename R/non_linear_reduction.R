#' Non linear dimensionality reduction for scRNAseq
#' 
#' t-SNE and UMAP functions are included in this package
#' 
#' @param input is the output of scaling (for function clusters_seurat()). Use the output of the clusters function as input for t-SNE and UMAP functions
#' @examples 
#' tsne_seurat(after_scaling.rds)
#' @export
# non-linear dimensionality reduction
clusters_seurat <- function(input){
   if(!(tools::file_ext(input)) == "rds") {
      return("Input file should be an rds file")
   } else if(tools::file_ext(input) == "rds") {
      clustering <- readRDS(input)
      
      if(class(clustering) != "Seurat") {
         return("File is not a seurat object")
      } else {
   clustering <- Seurat::RunPCA(clustering, features = Seurat::VariableFeatures(object = clustering))
   clustering <- Seurat::FindNeighbors(clustering, dims= 1:15)
   clustering <- Seurat::FindClusters(clustering, resolution = c(0.1, 0.3, 0.5, 0.7, 1))
   saveRDS(clustering, file = "clusters.rds")
      }
   }
   }

tsne_seurat <- function(input){
   if(!(tools::file_ext(input)) == "rds") {
      return("Input file should be an rds file")
   } else if(tools::file_ext(input) == "rds") {
      for_tsne <- readRDS(input)
      
      if(class(for_tsne) != "Seurat") {
         return("File is not a seurat object")
      } else {
   for_tsne <- Seurat::RunTSNE(for_tsne, dims = 1:10, dim.embed =2, label = TRUE)
   
   # load library for "element_text" to work properly
   library(ggplot2)
   library(plotly)
   
   new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
   names(new.cluster.ids) <- levels(for_tsne)
   for_tsne <- Seurat::RenameIdents(for_tsne, new.cluster.ids)
   tsne_plot <- Seurat::DimPlot(for_tsne, reduction = "tsne", label = TRUE) + xlab("t-SNE 1") + ylab("t-SNE 2") + theme(axis.title = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
   
   tsne_plot <- plotly::ggplotly(tsne_plot)
   
   htmltools::save_html(tsne_plot, file = "tsne_2dplot.html")
   
   
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
                               marker = list(size = 5, width=2), # controls size of points
                               text=~label, 
                               hoverinfo="text")
   
   htmltools::save_html(plotin3d, file = "tsne_3dplot.html")
   
      }
   }
   }

umap_seurat <- function(input){
   if(!(tools::file_ext(input)) == "rds") {
      return("Input file should be an rds file")
   } else if(tools::file_ext(input) == "rds") {
      for_umap <- readRDS(input)
      
      if(class(for_umap) != "Seurat") {
         return("File is not a seurat object")
      } else {
   for_3d <- for_umap
   
   # load library for "element_text" to work properly
   library(ggplot2)
   library(plotly)
   
   umap <- Seurat::RunUMAP(for_umap, dims = 1:10, n.components = 3L)
   new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
   names(new.cluster.ids) <- levels(umap)
   umap <- Seurat::RenameIdents(umap, new.cluster.ids)
   umap_plot <- Seurat::DimPlot(umap, reduction = "umap", label = TRUE) + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.title = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
   
   umap_plot <- plotly::ggplotly(umap_plot)
   htmltools::save_html(umap_plot, file = "umap_2dplot.html")
   
   # for 3d plot
   for_3d <- Seurat::RunUMAP(for_3d, dims = 1:10, n.components = 3L)
   plot.data <- Seurat::FetchData(object = for_3d, vars = c("umap_1", "umap_2", "umap_3", "seurat_clusters"))
   
   plot.data$label <- paste(rownames(plot.data))
   
   fig <- plot_ly(data = plot.data, 
                  x = ~umap_1, y = ~umap_2, z = ~umap_3, 
                  color = ~seurat_clusters,
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 5, width=2),
                  text=~label, 
                  hoverinfo="text")
   
   htmltools::save_html(fig, file = "umap_3dplot.html")
      }
   }
   }
