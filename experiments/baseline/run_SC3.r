suppressPackageStartupMessages({
  library(SC3)
  library(SingleCellExperiment)
  library(scater)
  library(aricode)
})


############################## Run sC3 ###################################

run_SC3 <- function(sce){
  
  if(ncol(sce) > 5000) {
    idx <- sample(1:ncol(sce), 5000)
    k <- sc3_estimate_k(sce[,idx])
    k <- as.numeric(metadata(k)$sc3)
    
    tmp <- sc3(sce, ks = k, gene_filter = F, svm_num_cells = 5000)
    tmp <- sc3_run_svm(tmp, ks = k)
  } else {
    tmp <- sc3(sce, k_estimator = T, gene_filter = F)
  }
  
  
  return(tmp)
}

#read data

############
#scDHA data#
############

#dataset <- 'tmuris'
#filepath <- paste0('~/R Scripts/rna_clustering/dataset/',dataset,'.rds')
#sce <- readRDS(filepath)
#label <- as.character(sce$cell_type1)

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

sce <- SingleCellExperiment(assays = list(logcounts = dt))
rowData(sce)$feature_symbol <- genes

label <- read.csv(labelspath)
label <- label[,1]


dt.SC3 <- run_SC3(sce)
  
cluster_class <- as.character(colData(dt.SC3)[,length(colData(dt.SC3))])
ARI(cluster_class, label)
NMI(cluster_class, label)

saveRDS(dt.SC3@metadata,'SC3_metadata_Mammary_Gland.rds')


