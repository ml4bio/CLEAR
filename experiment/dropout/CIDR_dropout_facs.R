suppressPackageStartupMessages({
  library(cidr)
  library(aricode)
  library(SingleCellExperiment)
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

#drop out rate
drop_rate <- 0.8

result_ARI <- c()
result_NMI <- c()

memory.limit(1e+10)

for (i in 1:5){
  sce_dropped <- dropout_sampling(dt, drop_rate = drop_rate, seed = i)
  sData <- run_CIDR(sce_dropped)
  
  ari <- ARI(sData@clusters, label)
  nmi <- NMI(sData@clusters, label)
  
  result_ARI <- c(ari, result_ARI)
  result_NMI <- c(nmi, result_NMI)
}

#mean & var
mean(result_ARI)
var(result_ARI)

mean(result_NMI)
var(result_NMI)
