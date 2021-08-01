# CLEAR

## Introduction

CLEAR (self-supervised **C**ontrastive **LEA**rning framework for sc**R**NA-seq) software provides a self-supervised learning based integrative single-cell RNA-seq data analysis tool. It overcomes the heterogeneity of the experimental data with a specifically designed representation learning task and thus can handle batch effects remove, dropout events, and time-trajectory inference better.

## Installation

The package can be installed by `git`. The test environment is CentOS 7 operation system, Nvidia TITAN X GPU.

### 1. Git clone from github

```bash
git clone https://github.com/ml4bio/CLEAR.git
cd ./CLEAR/
```

### 2. Use virtual environment with Anaconda
The main environment for CLEAR can be installed with this command:
```bash
conda env create -f environment.yml
conda activate CLEAR
```

## Running CLEAR

### 1. Prepare Data

To help users start their analysis quickly, we provide some python scripts, 
which can help users to convert the format of their files into h5ad format (required by CLEAR) and preprocess the data simultaneously.
For more individual preprocessing requirements, we also recommend some popular single-cell data analysis tools for data preprocessing instead of using the provided scripts.

#### (1). General Preprocessing

For those who need quick preparation for their data, we provide a Scanpy-based python script `generate_h5ad.py` in the `preprocess` folder to help with the format transformation and preprocessing.
It supports many kinds of input data formats, including 10X Genomics data, H5AD data, and CSV files. 
As for RDS files, you can transform them into the csv files by the script `rds_to_csv.py` in the `preprocess` folder and then preprocessed them with the same above script.
The required R environment of this R script needs to be set up by yourselves, of which the core package is `SingleCellExperiment`

If you only want to convert your processed data into the h5ad file, you can just set the input and output path when calling the script. For instance, when dealing with the h5ad file, you can run:
```bash
python preprocess/generate_h5ad.py --input_h5ad_path=Path_to_input --save_h5ad_dir=Path_to_Save_Folder
```
For 10X, you can use `--input_10X_path`; For csv, you can use `--count_csv_path` (and `--label_csv_path` for label if applicable)

Except for the format transformation, our script also provides many options for preprocessing, such as filtering, log, normalization, and so on. 
You can use `python generate_h5ad.py -h` for more details.

#### (2). Individual Preprocessing

For those who prefer a more individual data preparation, you can use [Scanpy](https://scanpy-tutorials.readthedocs.io/en/latest/) for preprocessing.
After completing all preprocessing steps, you can use the command `adata.write(USE_FOR_CLEAR.h5ad')` to save the h5ad file, and use this file as CLEAR input.

### 2. Apply CLEAR

With prepared input data, you can call the following script to implement CLEAR:
```bash
python CLEAR.py --input_h5ad_path="USE_FOR_CLEAR.h5ad" --epochs 100 --lr 0.01 --batch_size 512 --pcl_r 1024 --cos --gpu 0
```
Here, we only provide a set of commonly used CLEAR parameters for reference. You can run `python CLEAR.py -h` for more information.

**Note**: output files are saved in ./result/CLEAR, including `embeddings (feature.csv)`, `ground truth labels (if applicable)`, `cluster results (if applicable)` and some `log files (log)`.

You can then read the embeddings with Python (pd.read_csv) or R (read.csv) and incorperate it to the Anndata or Seurat for computing the neighborhood graph and following clustering.

## Running example

### 1. Download Dataset.

You can either download all of the datasets we used in our experiments with the script "download-data.sh" in the "data" folder, or use a single command to download specific dataset.

Here, we take "tabula-muris-senis-facs-processed-official-annotations-Bladder.h5ad " dataset for example.
```bash
wget https://ndownloader.figshare.com/files/23872610 -O data/original/h5ad/tmsfpoa-Bladder.h5ad
```

### 2. Generate Preprocessed H5AD File.
```bash
python preprocess/generate_h5ad.py --input_h5ad_path="./data/original/h5ad/tmsfpoa-Bladder.h5ad" --save_h5ad_dir="./data/preprocessed/h5ad/" --log
```

### 3. Apply CLEAR
```bash
python CLEAR.py --input_h5ad_path="./data/preprocessed/h5ad/tmsfpoa-Bladder_preprocessed.h5ad" --epochs 100 --lr 1 --batch_size 512 --pcl_r 1024 --cos --gpu 0
```
If you want to compute metric automatically with existing ground truth label, you can set `--obs_label_colname` to specify the column name of label.

## Citation

Han, W., Cheng, Y., Chen, J., Zhong, H., Hu, Z., Chen, S., . . . Li, Y. (2021). Self-supervised contrastive learning for integrative single cell RNA-seq data analysis. bioRxiv, 2021.2007.2026.453730. doi:10.1101/2021.07.26.453730

