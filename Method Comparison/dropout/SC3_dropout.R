######read files
suppressPackageStartupMessages({
  library(SC3)
  library(SingleCellExperiment)
  library(scater)
  library(aricode)
})

###########
#FACS data#
###########

#dataset <- 'Bladder'
#countspath <- paste0('~/R Scripts/rna_clustering/dataset/',dataset,'_counts.csv')
#labelspath <- paste0('~/R Scripts/rna_clustering/dataset/',dataset,'_celltype.csv')
#savepath <- paste0('~/R Scripts/rna_clustering/Seurat/Seurat_tabula-muris_',dataset,'.rds')
#save_label_path <- paste0('~/R Scripts/rna_clustering/labels_',dataset,'.rds')
#saveRDS(dt$cell_type1, save_label_path)
#dt <- read.csv(countspath)
#cells <- colnames(dt)[-1]
#genes <- dt[1][-1,]

#dt <- as.matrix(dt[-1,-1])

#row.names(dt) <- genes
#colnames(dt) <- cells

#label <- read.csv(labelspath)
#label <- label[,1]

#saveRDS(label, save_label_path)

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

############################# Run 5 trails ######################################

#read data

############
#scDHA data#
############

dataset <- 'lake'
filepath <- paste0('~/R Scripts/rna_clustering/dataset/',dataset,'.rds')
sce <- readRDS(filepath)
label <- as.character(sce$cell_type1)

#assay(sce,'logcounts') <- log10(assay(sce))

#drop out rate
#drop_rate <- 0.1

###########
#FACS data#
###########

dataset <- 'Diaphragm'
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

result_ARI <- c()
result_NMI <- c()

memory.limit(1e+10)

for (i in c(0.1,0.3,0.6,0.8)){
  drop_rate <- i
  sce_dropped <- dropout_sampling(sce, drop_rate = drop_rate, seed = 1)
  dt.SC3 <- run_SC3(sce_dropped)
  
  cluster_class <- as.character(colData(dt.SC3)[,length(colData(dt.SC3))])
  ari <- ARI(cluster_class, label)
  nmi <- NMI(cluster_class, label)
  
  name <- paste('~/R Scripts/rna_clustering/Dropout/SC3',i,dataset,'.rds')
  
  result_ARI <- c(result_ARI, ari)
  result_NMI <- c(result_NMI, nmi)
  
  saveRDS(dt.SC3@metadata, name)
}

result_ARI

#mean & var
mean(result_ARI)
var(result_ARI)

mean(result_NMI)
var(result_NMI)
