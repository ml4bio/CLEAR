suppressPackageStartupMessages({
  library(cidr)
  library(aricode)
  library(SingleCellExperiment)
})

memory.limit(1e+10)

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

############################## Run CIDR ###################################

run_CIDR <- function(dt){
  sData <- scDataConstructor(dt)
  sData <- determineDropoutCandidates(sData)
  sData <- wThreshold(sData)
  sData <- scDissim(sData)
  sData <- scPCA(sData)
  sData <- nPC(sData)
  
  sData <- scCluster(sData)
  
  return(sData)
}

############
#scDHA data#
############

dataset <- 'lake'
filepath <- paste0('~/R Scripts/rna_clustering/dataset/',dataset,'.rds')

sce <- readRDS(filepath)
cell_type1 <- as.character(sce$cell_type1)

#drop out rate
#drop_rate <- 0.8

###########
#FACS data#
###########

dataset <- 'Bladder'
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
  data <- assay(sce_dropped)
  sData <- run_CIDR(data)
  
  ari <- ARI(sData@clusters, label)
  nmi <- NMI(sData@clusters, label)
  
  result_ARI <- c(result_ARI, ari)
  result_NMI <- c(result_NMI, ari)
  
  name <- paste('~/R Scripts/rna_clustering/Dropout/CIDR',i,dataset,'.rds')
  saveRDS(sData,name)
}

result_ARI

#mean & var
mean(result_ARI)
var(result_ARI)

mean(result_NMI)
var(result_NMI)
