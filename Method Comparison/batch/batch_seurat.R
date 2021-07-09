######read files
suppressPackageStartupMessages({
  library(MatrixGenerics)
  library(Seurat)
  library(dplyr)
  library(SingleCellExperiment)
  library(aricode)
  library(mclust)
})

############################## Run seurat ###################################

run_Seurat <- function(sce){
  
  dt.seurat <- CreateSeuratObject(counts = sce, project = 'dataset_test')
  dt.seurat <- NormalizeData(dt.seurat, normalization.method = "LogNormalize", scale.factor = 10000, verbose = 0)
  dt.seurat <- FindVariableFeatures(dt.seurat, selection.method = "vst", nfeatures = 2000, verbose = 0)
  dt.seurat <- ScaleData(dt.seurat, verbose = 0)
  dt.seurat <- RunPCA(dt.seurat, features = VariableFeatures(object = dt.seurat), verbose = 0)
  dt.seurat <- FindNeighbors(dt.seurat, dims = 1:20, verbose = 0)
  dt.seurat <- FindClusters(dt.seurat, resolution = 0.5, verbose = 0)
  dt.seurat <- RunUMAP(dt.seurat, dims = 1:20, verbose = 0)
  
  return(dt.seurat)
}

#read data dataset 10
exp_matrix <- read.delim("rna_clustering/batch/dataset10/b1_exprs.txt", header = T, sep = "\t", )
exp_matrix[1:5,1:5]
rownames(exp_matrix) <- exp_matrix[,1]
exp_matrix <- exp_matrix[,-1]

exp_matrix_b2 <- read.delim("rna_clustering/batch/dataset10/b2_exprs.txt", header = T, sep = "\t", )
exp_matrix_b2[1:5,1:5]
rownames(exp_matrix_b2) <- exp_matrix_b2[,1]
exp_matrix_b2 <- exp_matrix_b2[,-1]

exp_matrix_total <- cbind(exp_matrix, exp_matrix_b2)

dt.seurat <- run_Seurat(exp_matrix_total)
dt.seurat$batch <- c(rep('b1',dim(exp_matrix)[2]), rep('b2', dim(exp_matrix_b2)[2]))

cluster_class <- as.character(dt.seurat@meta.data$seurat_clusters)

label_b1 <- read.delim('rna_clustering/batch/dataset10/b1_celltype.txt', header = T, sep = '\t')
label_b2 <- read.delim('rna_clustering/batch/dataset10/b2_celltype.txt', header = T, sep = '\t')

label_total <- rbind(label_b1, label_b2)
label_total$CellType
label_total$batch <- c(rep('b1',dim(exp_matrix)[2]), rep('b2', dim(exp_matrix_b2)[2]))

ARI(cluster_class, label_total$CellType)
NMI(cluster_class, label_total$CellType)

saveRDS(dt.seurat,'~/R Scripts/rna_clustering/Seurat/Seurat_dataset10.rds')
write.csv(cluster_class, 'predlabel_dataset10.csv')
write.csv(label_total, 'true_label.csv')

#read dataset 1
exp_matrix <- read.delim("rna_clustering/batch/dataset1/dataset1_sm_uc3.txt", header = T, sep = "\t", )
label <- read.delim('rna_clustering/batch/dataset1/sample_sm_uc3.txt', header = T, sep = '\t')
write.csv(label, 'true_label_dataset1.csv')

dt.seurat <- run_Seurat(exp_matrix)
dt.seurat$batch <- label$batch
cluster_class <- as.character(dt.seurat@meta.data$seurat_clusters)

ARI(cluster_class, label$celltype)
saveRDS(dt.seurat,'~/R Scripts/rna_clustering/Seurat/Seurat_dataset1.rds')

#read dataset 6
exp_matrix_b1 <- read.delim("rna_clustering/batch/dataset6/b1_exprs.txt", header = T, sep = "\t", )
exp_matrix_b1[1:5,1:5]
rownames(exp_matrix_b1) <- exp_matrix_b1[,1]
exp_matrix_b1 <- exp_matrix_b1[,-1]

exp_matrix_b2 <- read.delim("rna_clustering/batch/dataset6/b2_exprs.txt", header = T, sep = "\t", )
rownames(exp_matrix_b2) <- exp_matrix_b2[,1]
exp_matrix_b2 <- exp_matrix_b2[,-1]
exp_matrix_b2[1:5,1:5]

exp_matrix_b3 <- read.delim("rna_clustering/batch/dataset6/b3_exprs.txt", header = T, sep = "\t", )
rownames(exp_matrix_b3) <- exp_matrix_b3[,1]
exp_matrix_b3 <- exp_matrix_b3[,-1]
exp_matrix_b3[1:5,1:5]

exp_matrix_total <- cbind(exp_matrix_b1, exp_matrix_b2, exp_matrix_b3)

dt.seurat <- run_Seurat(exp_matrix_total)
dt.seurat$batch <- c(rep('b1',dim(exp_matrix_b1)[2]), 
                     rep('b2',dim(exp_matrix_b2)[2]), 
                     rep('b3',dim(exp_matrix_b3)[2]))

cluster_class <- as.character(dt.seurat$seurat_clusters)

label_b1 <- read.delim('rna_clustering/batch/dataset6/b1_celltype.txt', header = T, sep = '\t')
label_b2 <- read.delim('rna_clustering/batch/dataset6/b2_celltype.txt', header = T, sep = '\t')
label_b3 <- read.delim('rna_clustering/batch/dataset6/b3_celltype.txt', header = T, sep = '\t')

label_total <- rbind(label_b1, label_b2, label_b3)
label_total$CellType

saveRDS(dt.seurat,'~/R Scripts/rna_clustering/Seurat/Seurat_dataset6.rds')
write.csv(cluster_class, 'predlabel_dataset6.csv')
write.csv(label_total, 'true_label6.csv')

