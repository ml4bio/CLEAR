library(SingleCellExperiment)
library(scDHA)
library(aricode)

args<-commandArgs(TRUE)
print(args)
rds_root_path = args[1]
rds_name = args[2]
save_path = args[3]
print(rds_name)
#print(save_path)

# rds_root_path = "/home/yanhan/cjy/Single-Cell-Dataset/raw_rds/"
rds_path = paste(rds_root_path, rds_name, ".rds", sep = "")
# e.g., "/home/yanhan/cjy/Single-Cell-Dataset/raw_rds/baron-mouse.rds")
if(!file.exists(save_path)){
dir.create(save_path)
}
save_path = paste(save_path, "scDHA/", sep = "")
if(!file.exists(save_path)){
dir.create(save_path)
}
record_save_path = paste(save_path, "metric_scDHA.txt", sep = "")
latent_save_path = paste(save_path, "feature_scDHA_", rds_name ,".csv", sep = "")

file <- readRDS(rds_path)
data <- t(assay(file))
label <- as.character(file$cell_type1)

#Load example data (Goolam dataset)
#data("Goolam")

#Get data matrix and label
#data <- t(assay(Goolam))
#label <- as.character(Goolam$cell_type1)

#Log transform the data 
data <- log2(data + 1)

#Clustering
#Generate clustering result, the input matrix has rows as samples and columns as genes
result <- scDHA(data, seed = 1)

print(str(result))

#The clustering result can be found here 
cluster <- result$cluster

#Calculate adjusted Rand Index using mclust package
ari <- mclust::adjustedRandIndex(cluster,label)     # round(a, 2)
nmi <- NMI(cluster,label,variant = 'sum')
print(paste0("ARI: ", ari))
print(paste0("NMI: ", nmi))

#Save in
df = data.frame(
    DataSet = c(rds_name),
    ARI = c(ari),
    NMI = c(nmi)
)
write.table(df,file=record_save_path,append=T,row.names=F,col.names=F,quote=F)

#Save latent
latent <- result$latent
write.table(latent,file=latent_save_path,quote=F,row.names=F, col.names=F, sep=",")
#write.csv(latent,file=latent_save_path,quote=F,row.names = F) # can not omit col names
print("Successfully save data")

