######read files
suppressPackageStartupMessages({
  library(MatrixGenerics)
  library(Seurat)
  library(dplyr)
  library(SingleCellExperiment)
  library(aricode)
  library(mclust)
  library(scater)
})


################################ Drop out ####################################

#drop out function
dropout_sampling <- function(sce, drop_rate = 0, seed){
  #set seed
  set.seed(seed)
  #dropout sampling
  gene_size <- dim(sce)[1]
  drop <- sample(1:gene_size, as.integer(gene_size*(1-drop_rate)), replace = F)
  return(sce[drop,])
}

############################## Run seurat ###################################

run_Seurat <- function(sce){
  
  dt.seurat <- CreateSeuratObject(counts = assay(sce_dropped), project = 'dataset_test')
  dt.seurat <- NormalizeData(dt.seurat, normalization.method = "LogNormalize", scale.factor = 10000, verbose = 0)
  dt.seurat <- FindVariableFeatures(dt.seurat, selection.method = "vst", nfeatures = 2000, verbose = 0)
  dt.seurat <- ScaleData(dt.seurat, verbose = 0)
  dt.seurat <- RunPCA(dt.seurat, features = VariableFeatures(object = dt.seurat), verbose = 0)
  dt.seurat <- FindNeighbors(dt.seurat, dims = 1:20, verbose = 0)
  dt.seurat <- FindClusters(dt.seurat, resolution = 0.5, verbose = 0)
  dt.seurat <- RunUMAP(dt.seurat, dims = 1:20, verbose = 0)
  
  return(dt.seurat)
}

############################# Run 5 trails ######################################

#read data

############
#scDHA data#
############

dataset <- 'lake'
filepath <- paste0('~/R Scripts/rna_clustering/dataset/',dataset,'.rds')
savepath <- paste0('~/R Scripts/rna_clustering/Seurat/Seurat_',dataset,'.rds')
sce <- readRDS(filepath)
label <- as.character(sce$cell_type1)

###########
#FACS data#
###########

dataset <- 'Mammary_Gland'
countspath <- paste0('~/R Scripts/rna_clustering/dataset/',dataset,'_counts.csv')
labelspath <- paste0('~/R Scripts/rna_clustering/dataset/',dataset,'_celltype.csv')

dt <- read.csv(countspath)
cells <- colnames(dt)[-1]
genes <- dt[1][-1,]

dt <- as.matrix(dt[-1,-1])

rownames(dt) <- genes
colnames(dt) <- cells

sce <- SingleCellExperiment(assays = list(counts = dt))
sce <- logNormCounts(sce)

rowData(sce)$feature_symbol <- genes

label <- read.csv(labelspath)
label <- label[,1]

#drop out rate
#drop_rate <- 0.1

result_ARI <- c()
result_NMI <- c()

memory.limit(1e+10)

for (i in c(0.1,0.3,0.6,0.8)){
  drop_rate <- i
  sce_dropped <- dropout_sampling(sce, drop_rate = drop_rate, seed = 1)
  dt.seurat <- run_Seurat(sce_dropped)
  
  cluster_class <- as.character(dt.seurat@meta.data$seurat_clusters)
  ari <- ARI(cluster_class, label)
  nmi <- NMI(cluster_class, label)
  
  result_ARI <- c(result_ARI, ari)
  result_NMI <- c(result_NMI, nmi)
  
  name <- paste('~/R Scripts/rna_clustering/Dropout/Seurat',i,dataset,'.rds')
  
  saveRDS(dt.seurat@reductions,name)
}

result_ARI

#mean & var
mean(result_ARI)
var(result_ARI)

mean(result_NMI)
var(result_NMI)
