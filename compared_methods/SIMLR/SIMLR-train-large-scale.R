# required external packages for SIMLR large scale
library("Rcpp")
library("Matrix")
library("pracma")
library("RcppAnnoy")
library("RSpectra")


# load the igraph package to compute the NMI
library("igraph")

# load the palettes for the plots
library(grDevices)

library(aricode)

# load the SIMLR R package
source('./method/SIMLR/R/SIMLR_Large_Scale.R')
source("./method/SIMLR/R/utils.simlr.large.scale.R")
source("./method/SIMLR/R/utils.simlr.R")
source("./method/SIMLR/R/SIMLR.Rtsne.R")

# load the C file

# NOTE 1: we make use of an external C program during the computations of SIMLR.
# The code is located in the R directory in the file projsplx_R.c. In order to 
# use SIMLR one needs to compite the program. To do so, one needs to run on the 
# shell the command R CMD SHLIB -c projsplx_R.c. 
# The precompiled projsplx_R.so is already provided for MAC OS X only. 
# If one wants to use SIMLR on other operative systems, the file projsplx_R.so 
# needs to be deleted, and re-compiled. 

# NOTE 2: for Windows, the command dyn.load("./R/projsplx_R.so") needs to be 
# substituted with the command dyn.load("./R/projsplx_R.dll"). 

dyn.load("./method/SIMLR/R/projsplx_R.so")

print("A")

sourceCpp("./method/SIMLR/src/Rtsne.cpp")



args<-commandArgs(TRUE)
print(args)
csv_root_path = args[1]
dataset_name = args[2]
save_path = args[3]
print(dataset_name)
#print(save_path)

# csv_root_path = "/home/yanhan/cjy/Single-Cell-Dataset/raw_rds/"
csv_count_path = paste(csv_root_path, dataset_name, "_counts.csv", sep = "")
csv_label_path = paste(csv_root_path, dataset_name, "_labels.csv", sep = "")
# e.g., "/home/yanhan/cjy/Single-Cell-Dataset/raw_rds/baron-mouse.rds")
if(!file.exists(save_path)){
dir.create(save_path)
}
save_path = paste(save_path, "SIMLR/", sep = "")
if(!file.exists(save_path)){
dir.create(save_path)
}
record_save_path = paste(save_path, "metric_SIMLR.txt", sep = "")
latent_save_path = paste(save_path, "feature_SIMLR_", dataset_name ,".csv", sep = "")
gt_label_save_path = paste(save_path, "gt_label_SIMLR_", dataset_name ,".csv", sep = "")
pd_label_save_path = paste(save_path, "pd_label_SIMLR_", dataset_name ,".csv", sep = "")
pred_save_path = paste(save_path, "pred_SIMLR_", dataset_name ,".csv", sep = "")

# discard for rds
#file <- readRDS(rds_path)
#data <- t(assay(file))
#label <- as.character(file$cell_type1)

# rewrite for csv
dt <- read.csv(csv_count_path)  # row - cell; col - gene
cells <- dt[1][-1,]
genes <- colnames(dt)[-1]
#cells <- colnames(dt)[-1]
#genes <- dt[1][-1,]
#print(str(dt))
dt <- as.matrix(dt[,-1])
#row.names(dt) <- genes
#colnames(dt) <- cells
data <- dt
data = t(data)
label <- read.csv(csv_label_path)[,2]
print(dim(data))
print(length(label))
print(str(data))

#Clustering
clust_list = unique(label)
n_clust = length(clust_list)
print(clust_list)
print(n_clust)

set.seed(11111)
result = SIMLR_Large_Scale(X=data,c=n_clust)#,k=30,kk=200)
#print(str(result))

#The clustering result can be found here
cluster <- result$y$cluster

print("cluster")
print(str(cluster))
print("label")
print(str(label))

#Calculate adjusted Rand Index using mclust package
ari <- mclust::adjustedRandIndex(cluster,label)     # round(a, 2)
nmi <- NMI(cluster,label,variant = 'sum')
print(paste0("ARI: ", ari))
print(paste0("NMI: ", nmi))

#Save in
df = data.frame(
    DataSet = c(dataset_name),
    ARI = c(ari),
    NMI = c(nmi)
)
write.table(df,file=record_save_path,append=T,row.names=F,col.names=F,quote=F)

#Save latent
latent <- result$latent
write.table(latent,file=latent_save_path,quote=F,row.names=F, col.names=F, sep=",")
#write.csv(latent,file=latent_save_path,quote=F,row.names = F) # can not omit col names
print("Successfully save data")

# Save label
write.csv(label,file=gt_label_save_path, quote=F,row.names = T)

write.csv(cluster,file=pd_label_save_path, quote=F,row.names = T)



