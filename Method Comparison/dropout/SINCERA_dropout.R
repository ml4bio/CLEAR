suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SINCERA)
  library(aricode)
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

################################ Run SINCERA ####################################

run_SINCERA <- function(data){
  ###Construct S4 object for sincera
  data <- as.data.frame(data)
  sc <- construct(exprmatrix = data, samplevector = paste("sample", rep(1,length(cell_type1)), sep=""))
  
  ###Apply hierarchy cluster
  
  sc <- normalization.zscore(sc, pergroup=TRUE)
  
  sc <- cluster.assignment(sc,clustering.method="hc", distance.method="pearson", verbose = F)
  
  cluster <- sc@data@phenoData@data$CLUSTER
  
  return(cluster)
}

############
#scDHA data#
############

dataset <- 'lake'
filepath <- paste0('~/R Scripts/rna_clustering/dataset/',dataset,'.rds')

sce <- readRDS(filepath)
cell_type1 <- as.character(sce$cell_type1)

#drop out rate
drop_rate <- 0.8

result_ARI <- c()
result_NMI <- c()

memory.limit(1e+10)

for (i in 1:5){
  sce_dropped <- dropout_sampling(sce, drop_rate = drop_rate, seed = i)
  data <- assay(sce_dropped)
  cluster <- run_SINCERA(data)
  
  ari <- ARI(cluster, cell_type1)
  nmi <- NMI(cluster, cell_type1)
  
  result_ARI <- c(ari, result_ARI)
  result_NMI <- c(nmi, result_NMI)
}

#mean & var
mean(result_ARI)
var(result_ARI)

mean(result_NMI)
var(result_NMI)