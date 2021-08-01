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

### 1. Prepare Dataset

To help users start their analysis quickly, we provide some python scripts to use. For more individual preprocessing requirements, we also recommend some popular single-cell data analysis tools for data preprocessing instead of using the provided scripts.

#### (1). 10X Genomics data

The widely-used 10X platform provides a `10x matrix` data structure: `barcodes.tsv.gz, features.tsv.gz, matrix.tsv.gz`

i. (Recommanded) Use [Scanpy](https://scanpy-tutorials.readthedocs.io/) for data loading and individual data filtering.

```python
import scanpy as sc
adata = sc.read_10x_mtx(data_path)
...(data preprocessing steps with individual paramters, see Scanpy tutorial)
adata.write('USE_FOR_CLEAR.h5ad')
```

ii. Use `10X_to_h5ad.py' for quick start. `10X_to_h5ad.py' provides a more uniform data filtering function.

```python

```

#### (2). CSV file

If a expression profile (usually a csv file) is provided, 

i. (Recommanded) Construct `Anndata` and preprocess the data by yourself with more unique filtering parameters

```python
#create anndata
df = pd.read_csv('count.csv', index_col=0)
cellinfo = pd.DataFrame(df.index,index=df.index,columns=['sample_index'])
geneinfo = pd.DataFrame(df.columns,index=df.columns,columns=['genes_index'])
adata = sc.AnnData(df, obs=cellinfo, var = geneinfo)

#generate h5ad file for CLEAR input
...(data preprocessing steps with individual paramters, see Scanpy tutorial)
adata.write('USE_FOR_CLEAR.h5ad')
```

ii. Use `csv_to_h5ad.py` for quick start. `csv_to_h5ad.py` provides a more uniform data filtering function. 

```python

```

#### (3). Rdata

For R tool users, we highly recommand you save the filtered csv file and create a h5ad file as input.

```python
#create anndata
df = pd.read_csv('R_FILTERED_DATA.csv', index_col=0)
cellinfo = pd.DataFrame(df.index,index=df.index,columns=['sample_index'])
geneinfo = pd.DataFrame(df.columns,index=df.columns,columns=['genes_index'])
adata = sc.AnnData(df, obs=cellinfo, var = geneinfo)

#generate h5ad file for CLEAR input
adata.write('USE_FOR_CLEAR.h5ad')
```


### 2. Apply CLEAR

we can apply CLEAR with the following command:
```
python CLEAR.py --input_h5ad_path="USE_FOR_CLEAR.h5ad" --obs_label_colname="cell_ontology_class" --epochs 100 --lr 0.01 --batch_size 512 --pcl_r 1024 --cos --gpu 0
```
Note: output files are saved in ./result/CLEAR, including embeddings (feature.csv), ground truth labels (gt_label.csv if applicable), cluster results (pd_label.csv if applicable) and some log files (log)


## Citation

Han, W., Cheng, Y., Chen, J., Zhong, H., Hu, Z., Chen, S., . . . Li, Y. (2021). Self-supervised contrastive learning for integrative single cell RNA-seq data analysis. bioRxiv, 2021.2007.2026.453730. doi:10.1101/2021.07.26.453730

