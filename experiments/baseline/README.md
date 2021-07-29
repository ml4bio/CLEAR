# Benchmark Experiments
There are procedures for our method, CLEAR, and other compared methods implemented on the benchmark datasets.

Some of the methods can be implemented under the following instructions
## 1.Preprocess
First, we can preporcess the datasets with our scripts in the "preprocess" folder.

## 2.Run Methods
CLEAR (ours)
```
python method/Scarlet/main_pcl_cjy.py --lr 10 --batch-size 512 --pcl-r 1024 --cos --count_data "${csv_dir}${dataset_name}_counts.csv" --label_data "${csv_dir}${dataset_name}_labels.csv" --exp-dir exp_tmp --gpu 0 --epochs 500
```
scVI, scDHA, ItClust, scGNN, SIMLR
```
python method/scVI-train.py ${h5ad_dir} ${dataset_name} ${save_dir}
Rscript method/scDHA-train.R ${csv_dir} ${dataset_name} ${save_dir}
python method/ItClust-train.py ${h5ad_dir} ${dataset_name} ${save_dir}
python method/scGNN/scGNN-cluster.py ${csv_dir} ${dataset_name} ${save_dir}
Rscript method/SIMLR/SIMLR-train.R ${csv_dir} ${dataset_name} ${save_dir}
```

while CIDR,SC3,Seurat,SINCERA should be directly applied with the corresponding R scripts, which already include the preprocessing step.