# non-linear-reduction
This package contains three functions: a) `clusters_seurat()`, b) `tsne_seurat()`, and c) `umap_seurat`.

Seurat utilizes a graph-based clustering approach by applying modularity optimization techniques to iteratively group cells together. This is naturally implemented by a Seurat function deployed by the `clusters_seurat()` function.

The last two functions, `tsne_seurat()` and `umap_seurat()` deploy non-linear dimensional reduction techniques to visualize and explore scRNAseq data. 

## Installation
The package can be installed using
```
devtools::install_github("BTIP/non-linear-reduction")
```

## Example
The output of the `clusters_seurat()` is an rds file that functions as an input file for the remaining two functions. These last two functions would produce 2D and 3D plots.
```
# to perform clustering
clusters_seurat("after_scaling.rds")

# use the same output to the two functions
tsne_seurat("clusters.rds")

umap_seurat("clusters.rds")
```
## scRNAseq processing workflow 

The standard scRNAseq processing workflow with the R package Seurat consists of seven (7) steps. This package is the last package that contains the remaining steps of the pipeline.

The following are the repositories of the packages for every step of the pipeline:
1. QC and filtering: [qualitycontrolseurat package](https://github.com/BTIP2024/quality-control-seurat)
2. Normalization: [qualitycontrolseurat package](https://github.com/BTIP2024/quality-control-seurat)
3. Identification of highly variable features: [selectionscalingseurat package](https://github.com/BTIP2024/selection-scaling-seurat)
4. Scaling: [selectionscalingseurat package](https://github.com/BTIP2024/selection-scaling-seurat)
5. Linear Dimensionality Reduction (PCA): [pcaseurat package](https://github.com/BTIP2024/pca-seurat)
6. Clustering: [nonlinearreduction package](https://github.com/BTIP2024/non-linear-reduction)
7. Non-linear dimensionality reduction (t-SNE and UMAP): [nonlinearreduction package](https://github.com/BTIP2024/non-linear-reduction)

An overview of the pipeline and its outputs can be observed below:
![](https://github.com/user-attachments/assets/fed7196a-bce5-447e-96e7-eaaef392e8d6)
