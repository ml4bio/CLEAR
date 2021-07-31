# CLEAR

## Introduction

CLEAR (self-supervised **C**ontrastive **LEA**rning framework for sc**R**NA-seq) software provides a self-supervised learning based integrative single-cell RNA-seq data analysis tool. It overcomes the heterogeneity of the experimental data with a specifically designed representation learning task and thus can handle batch effects remove, dropout events, and time-trajectory inference better.

## Installation

The package can be installed based by `git `. Test environment is CentOS 7 operation system, Nvidia TITAN X GPU.

### 1. (Recommanded) Use virtual environment with Anaconda
The main environment for CLEAR can be installed with this command:
```
conda env create -f CLEAR_environment.yml

```
To process rds data, we should also create a R environment:
```
conda env create -f Rdata_environment.yml
```

### 2. Git clone from github

```
git clone https://github.com/ml4bio/CLEAR.git
cd ~/CLEAR/
```

## Quick Running

### 1. Prepare Dataset

There are two kinds of input data format: rds and h5ad. The preprocessing step will , so we should transform them into csv files respectively. 
In the following examples, I will use baron-mouse.rds and abula-muris-senis-facs-processed-official-annotations-Diaphragm.h5ad as references. 
You can either download all of them with the script "download-data.sh" in the "data" folder or use the command in it to download specific dataset.
Here, we take "deng.rds" dataset for example.

(1) download dataset.
```
wget https://bioinformatics.cse.unr.edu/software/scDHA/resource/Reproducibility/Data/deng.rds -O data/original/rds/deng.rds
wget https://ndownloader.figshare.com/files/23872487 -O data/original/h5ad/tmsfpoa-Diaphragm.h5ad
```

(2) convert format, from rds files to csv files.

We can use the offered script "rds_to_csv.R" in the "preprocess" folder to accomplish this transfomation. The command line is as follows: 
```
conda activate Rdata
Rscript preprocess/rds_to_csv.R ./data/original/rds/ deng data/original/csv/
```
Then you can find the baron-mouse_counts.csv and baron-mouse_labels.csv inside the data/original/csv folder.
As for input h5ad files, you can use preprocess/h5ad_to_csv.py to convert them into csv format.
```
conda activate CLEAR
python preprocess/h5ad_to_csv.py "./data/original/h5ad/abula-muris-senis-facs-processed-official-annotations-Diaphragm.h5ad" "./data/original/csv/" 0 0 0 0
```

(3) preprocess csv files and generate input h5ad file.

```
conda activate CLEAR
python preprocess/preprocess_csv_to_h5ad.py --count_csv_path="./data/original/csv/deng_counts.csv" --label_csv_path="./data/original/csv/deng_labels.csv" --save_h5ad_dir="./data/preprocessed/h5ad/" --label_colname="x" --log --drop_prob=0

python preprocess/generate_preprocessed_h5ad.py --input_h5ad_path="./data/original/h5ad/tmsfpoa-Diaphragm.h5ad" --save_h5ad_dir="./data/preprocessed/h5ad/" --label_colname="x" --log --drop_prob=0
```

### 2. Apply CLEAR

we can apply CLEAR with the following command:
```
python clear/clear.py --input_h5ad_path="./data/preprocessed/h5ad/deng.h5ad" --epochs 100 --lr 1 --batch_size 512 --pcl_r 1024 --cos --gpu 0

python CLEAR.py --input_h5ad_path="./data/preprocessed/h5ad/tmsfpoa-Diaphragm_preprocessed.h5ad" --epochs 100 --lr 0.01 --batch_size 512 --pcl_r 1024 --cos --gpu 0
```
Note: output files are saved in ./result/CLEAR, including embeddings, ground truth labels, cluster results and some log files


## Citation

Han, W., Cheng, Y., Chen, J., Zhong, H., Hu, Z., Chen, S., . . . Li, Y. (2021). Self-supervised contrastive learning for integrative single cell RNA-seq data analysis. bioRxiv, 2021.2007.2026.453730. doi:10.1101/2021.07.26.453730

