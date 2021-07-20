library(SingleCellExperiment)

print("rds to csv")

args<-commandArgs(TRUE)
#print(args)

root_path = args[1]
rds_name = args[2]
save_path = args[3]
print(rds_name)
#print(save_path)

rds_path = paste(root_path, rds_name, ".rds", sep = "")
if(!file.exists(rds_path)){
rds_path = paste(root_path, rds_name, ".RDS", sep = "")
}
# e.g., "/home/yanhan/cjy/Single-Cell-Dataset/raw_rds/baron-mouse.rds")

file <- readRDS(rds_path)
data <- t(assay(file))
label <- as.character(file$cell_type1)

# if you want to save the file as csv
if(!file.exists(save_path)){
dir.create(save_path)
}

save_data_path = paste(save_path, rds_name, "_counts.csv", sep = "")
save_label_path = paste(save_path, rds_name, "_labels.csv", sep = "")

print(save_data_path)
print(save_label_path)
write.csv(data, save_data_path)
write.csv(label, save_label_path)
print("complete")
