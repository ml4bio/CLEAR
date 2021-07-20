# for dataset1, dataset6 apply other method
h5ad_dir="/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/h5ad/"
csv_dir="/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/csv/"
save_dir="/home/yanhan/cjy/Single-Cell-Dataset/Single-Cell-Cluster/result_batch_effect/"

dataset_name="dataset1"
#python method/scVI-train.py ${h5ad_dir} ${dataset_name} ${save_dir}
#Rscript method/scDHA-train.R ${csv_dir} ${dataset_name} ${save_dir}
#python method/ItClust-train.py ${h5ad_dir} ${dataset_name} ${save_dir}
Rscript method/SIMLR/SIMLR-train.R ${csv_dir} ${dataset_name} ${save_dir}

dataset_name="dataset6"
#python method/scVI-train.py ${h5ad_dir} ${dataset_name} ${save_dir}
#Rscript method/scDHA-train.R ${csv_dir} ${dataset_name} ${save_dir}
#python method/ItClust-train.py ${h5ad_dir} ${dataset_name} ${save_dir}
Rscript method/SIMLR/SIMLR-train.R ${csv_dir} ${dataset_name} ${save_dir}


# dataset3
h5ad_dir="/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/h5ad/"
csv_dir="/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/csv/"
save_dir="/home/yanhan/cjy/Single-Cell-Dataset/Single-Cell-Cluster/result_batch_effect/"
dataset_names=("dataset3_simul1" "dataset3_simul1_HVG" "dataset3_simul2" "dataset3_simul2_HVG" "dataset3_simul3" "dataset3_simul3_HVG" "dataset3_simul4" "dataset3_simul4_HVG" "dataset3_simul5" "dataset3_simul5_HVG" "dataset3_simul6" "dataset3_simul6_HVG")

for dataset_name in ${dataset_names[@]};
do
  #python method/scVI-train.py ${h5ad_dir} ${dataset_name} ${save_dir}
  #Rscript method/scDHA-train.R ${csv_dir} ${dataset_name} ${save_dir}
  #python method/ItClust-train.py ${h5ad_dir} ${dataset_name} ${save_dir}
  Rscript method/SIMLR/SIMLR-train.R ${csv_dir} ${dataset_name} ${save_dir}
done