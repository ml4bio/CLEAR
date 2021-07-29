suppressPackageStartupMessages({
  library(cidr)
  library(aricode)
  library(SingleCellExperiment)
})

set.seed(2021)

memory.limit(1e11)

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

#dataset <- 'lake'
#filepath <- paste0('~/R Scripts/rna_clustering/dataset/',dataset,'.rds')

#sce <- readRDS(filepath)
#cell_type1 <- as.character(sce$cell_type1)

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

row.names(dt) <- genes
colnames(dt) <- cells

label <- read.csv(labelspath)[,1]

sData <- run_CIDR(dt)
  
ARI(sData@clusters, label)
NMI(sData@clusters, label)

saveRDS(sData, '~/R Scripts/rna_clustering/CIDR/CIDR_sData_Mammary_Gland.rds')
