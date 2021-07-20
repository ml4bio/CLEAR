######read files
suppressPackageStartupMessages({
  library(MatrixGenerics)
  library(Seurat)
  library(dplyr)
  library(SingleCellExperiment)
  library(aricode)
  library(mclust)
})

memory.limit(1e5)

############################## Run seurat ###################################

run_Seurat <- function(sce){
  
  dt.seurat <- CreateSeuratObject(counts = sce, project = 'dataset_test', min.cells = 3, min.features = 200)
  dt.seurat <- NormalizeData(dt.seurat, normalization.method = "LogNormalize", scale.factor = 10000, verbose = 0)
  dt.seurat <- FindVariableFeatures(dt.seurat, selection.method = "vst", nfeatures = 2000, verbose = 0)
  dt.seurat <- ScaleData(dt.seurat, verbose = 0)
  dt.seurat <- RunPCA(dt.seurat, features = VariableFeatures(object = dt.seurat), verbose = 0)
  dt.seurat <- FindNeighbors(dt.seurat, dims = 1:20, verbose = 0)
  dt.seurat <- FindClusters(dt.seurat, resolution = 0.5, verbose = 0)
  dt.seurat <- RunUMAP(dt.seurat, dims = 1:20, verbose = 0)
  
  return(dt.seurat)
}

#read data

############
#scDHA data#
############

#dataset <- 'pollen'
#filepath <- paste0('~/R Scripts/rna_clustering/dataset/',dataset,'.rds')
savepath <- paste0('~/R Scripts/rna_clustering/Seurat/Seurat_',dataset,'.csv')
#sce <- readRDS(filepath)
#label <- as.character(sce$cell_type1)

labelepath <- paste0('~/R Scripts/rna_clustering/Seurat/Seurat_label_',dataset,'.csv')
write.csv(label, labelepath)
###########
#FACS data#
###########

dataset <- 'Limb_Muscle'
countspath <- paste0('~/R Scripts/rna_clustering/dataset/',dataset,'_counts.csv')
labelspath <- paste0('~/R Scripts/rna_clustering/dataset/',dataset,'_celltype.csv')

dt <- read.csv(countspath)
cells <- colnames(dt)[-1]
genes <- dt[1][-1,]

dt <- as.matrix(dt[-1,-1])

row.names(dt) <- genes
colnames(dt) <- cells

label <- read.csv(labelspath)
label <- label[,1]

dt.seurat <- run_Seurat(dt)
cluster_class <- as.character(dt.seurat@meta.data$seurat_clusters)

#dt.seurat@reductions

ARI(cluster_class, label)
NMI(cluster_class, label)

write.csv(cluster_class, savepath)
saveRDS(dt.seurat@reductions,'~/R Scripts/rna_clustering/Seurat/Seurat_reducation_Mammary.rds')
