library(SingleCellExperiment)
library(scDHA)
library(aricode)

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
save_path = paste(save_path, "scDHA/", sep = "")
if(!file.exists(save_path)){
dir.create(save_path)
}
record_save_path = paste(save_path, "metric_scDHA.txt", sep = "")
latent_save_path = paste(save_path, "feature_scDHA_", dataset_name ,".csv", sep = "")
gt_label_save_path = paste(save_path, "gt_label_scDHA_", dataset_name ,".csv", sep = "")
pd_label_save_path = paste(save_path, "pd_label_scDHA_", dataset_name ,".csv", sep = "")
pred_save_path = paste(save_path, "pred_scDHA_", dataset_name ,".csv", sep = "")

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
label <- read.csv(csv_label_path)[,2]
print(dim(data))
print(length(label))

#Log transform the data.
data <- log2(data + 1)

#Clustering
#Generate clustering result, the input matrix has rows as samples and columns as genes
result <- scDHA(data, seed = 25)  #1

# 2d-Visualization
#Generate 2D representation, the input is the output from scDHA function
result <- scDHA.vis(result, seed = 1)
#Plot the representation of the dataset, different colors represent different cell types
#plot(result$pred, col=factor(label), xlab = "scDHA1", ylab = "scDHA2")
# save
write.table(result$pred,file=pred_save_path,quote=F,row.names=F, col.names=F, sep=",")

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

