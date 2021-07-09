suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SINCERA)
  library(aricode)
})

set.seed(2021)

memory.limit(1e11)

################################ Run SINCERA ####################################

run_SINCERA <- function(data){
  ###Construct S4 object for sincera
  data <- as.data.frame(data)
  sc <- construct(exprmatrix = data, samplevector = paste("sample", rep(1,length(cell_type1)), sep=""))
  
  ###Apply hierarchy cluster
  
  sc <- normalization.zscore(sc, pergroup=TRUE)
  
  sc <- cluster.assignment(sc,clustering.method="hc", distance.method="pearson", verbose = F)
  
  return(sc)
}

############
#scDHA data#
############

dataset <- 'yan'
filepath <- paste0('~/R Scripts/rna_clustering/dataset/',dataset,'.rds')

sce <- readRDS(filepath)
cell_type1 <- as.character(sce$cell_type1)


data <- assay(sce)
sc <- run_SINCERA(data)

cluster <- sc@data@phenoData@data$CLUSTER

ARI(cluster, cell_type1)
NMI(cluster, cell_type1)


saveRDS(cluster, 'SINCERA_sc_yan.rds')
