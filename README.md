# CLEAR

## Introduction

CLEAR (self-supervised **C**ontrastive **LEA**rning framework for sc**R**NA-seq) software provides a self-supervised learning based integrative single-cell RNA-seq data analysis tool. It overcomes the heterogeneity of the experimental data with a specifically designed representation learning task and thus can handle batch effects remove, dropout events, and time-trajectory inference better.

## Installation

The package can be installed based by `git `. Test environment is CentOS 7 operation system, Nvidia TITAN X GPU.

### 1. Git clone from github

```
git clone https://github.com/ml4bio/CLEAR.git
cd ~/CLEAR/
```

### 2. Use virtual environment with Anaconda
The main environment for CLEAR can be installed with this command:
```
conda env create -f environment.yml
```

## Running CLEAR

### 1. Prepare Data

To help users start their analysis quickly, we provide some python scripts to use. 
For more individual preprocessing requirements, we also recommend some popular single-cell data analysis tools for data preprocessing instead of using the provided scripts.

#### (1). Only Convert Format

For H5AD data,
```bash
python preprocess/generate_h5ad.py --input_h5ad_path=Path_to_10X --save_h5ad_dir=Path_to_Save_Folder
```

For 10X Genomics data

The widely-used 10X platform provides a `10x matrix` data structure: `barcodes.tsv.gz, features.tsv.gz, matrix.tsv.gz`
```bash
python preprocess/generate_h5ad.py --input_10X_path=Path_to_10X --save_h5ad_dir=Path_to_Save_Folder
```

For CSV data,  If a expression profile (usually a csv file) is provided, 

```bash
python preprocess/generate_h5ad.py --count_csv_path=Path_to_Count_CSV --label_csv_path=Path_to_Label_CSV --save_h5ad_dir=Path_to_Save_Folder
```

As for rds files, you can transform them into csv files by the script "rds_to_csv.py" in the "preprocess" folder and then preprocessed them with the same above script.
The required R environment of this R script need to be set up by yourselves, of which the core package is "SingleCellExperiment"

#### (2). Preprocess if needed
You can use your own script to accompulish preprocessing and you can also add some input parameters when using the script we offered above. 



### 2. Apply CLEAR

we can apply CLEAR with the following command:
```
python CLEAR.py --input_h5ad_path="USE_FOR_CLEAR.h5ad" --epochs 100 --lr 0.01 --batch_size 512 --pcl_r 1024 --cos --gpu 0
```
Note: output files are saved in ./result/CLEAR, including embeddings (feature.csv), ground truth labels (gt_label.csv if applicable), cluster results (pd_label.csv if applicable) and some log files (log)

## Example
Here, we take "deng.rds" dataset for example.
### 1.download dataset.
You can either download all of them with the script "download-data.sh" in the "data" folder or use the command in it to download specific dataset.
Here, we take "deng.rds" dataset for example.
```
wget https://ndownloader.figshare.com/files/23872610 -O data/original/h5ad/tmsfpoa-Bladder.h5ad
```

### 2.generate preprocessed h5ad file.
```
python preprocess/generate_h5ad.py --input_h5ad_path="./data/original/h5ad/tmsfpoa-Bladder.h5ad" --save_h5ad_dir="./data/preprocessed/h5ad/" --log
```

### 3.apply CLEAR
```
python CLEAR.py --input_h5ad_path="./data/preprocessed/h5ad/tmsfpoa-Bladder_preprocessed.h5ad" --obs_label_colname="cell_ontology_class" --epochs 100 --lr 1 --batch_size 512 --pcl_r 1024 --cos --gpu 0
```

## Citation

Han, W., Cheng, Y., Chen, J., Zhong, H., Hu, Z., Chen, S., . . . Li, Y. (2021). Self-supervised contrastive learning for integrative single cell RNA-seq data analysis. bioRxiv, 2021.2007.2026.453730. doi:10.1101/2021.07.26.453730

