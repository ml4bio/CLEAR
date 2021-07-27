# CLEAR

## Introduction

CLEAR (self-supervised **C**ontrastive **LEA**rning framework for sc**R**NA-seq) software provides a self-supervised learning based integrative single-cell RNA-seq data analysis tool. It overcomes the heterogeneity of the experimental data with a specifically designed representation learning task and thus can handle batch effects remove, dropout events, and time-trajectory inference better.

## Installation

The package can be installed based by `git `. Test environment is CentOS 7 operation system, Nvidia TITAN X GPU.

### (Recommanded) Use virtual environment with Anaconda (Python = 3.7)

```bash
conda create -n CLEAR python = 3.7
conda activate CLEAR
```

### Git from github

```bash
git https://github.com/ml4bio/CLEAR
cd ~/CLEAR/CLEAR
```

## Quick Running

1.preprocess datasets

There are two kinds of input data format: rds and h5ad. The required of our script is a pair of csv files including "counts.csv" and "labels.csv"., so we should transform them into csv files respectively. In the following examples, I will use baron-mouse.rds and abula-muris-senis-facs-processed-official-annotations-Diaphragm.h5ad as references. You can download them with the script "download-data.sh" in the "data" folder.

(1) rds files to csv files

We can use the offered script "rds_to_csv.R" in the "preprocess" folder to accomplish this transfomation. The command line is as follows: 
```
Rscript preprocess/rds_to_csv.R data/rds/baron-mouse.rds data/ocsv/
```
Then you can find the baron-mouse_counts.csv and baron-mouse_labels.csv inside the data/ocsv folder.

(2) h5ad files to csv files (without preprocess)
```
python preprocess/h5ad_to_csv.py data/oh5ad/abula-muris-senis-facs-processed-official-annotations-Diaphragm.h5ad data/ocsv 0 0
```


2. apply CLEAR

we can apply CLEAR with the following command:

For baron-mouse,
```
python /CLEAR/main_pcl_cjy.py --lr 10 --batch-size 512 --pcl-r 1024 --cos --count_data "data/ocsv/baron-mouse_counts.csv" --label_data "data/ocsv/baron-mouse_labels.csv" --exp-dir CLEAR/exp_tmp --gpu 0 --epochs 500
```

For abula-muris-senis-facs-processed-official-annotations-Diaphragm
```
python method/Scarlet/main_pcl_cjy.py --lr 5e-1 --batch-size 512 --pcl-r 1024 --cos --count_data "data/ocsv/abula-muris-senis-facs-processed-official-annotations-Diaphragm_counts.csv" --label_data "data/ocsv/abula-muris-senis-facs-processed-official-annotations-Diaphragm_labels.csv" --exp-dir exp_tmp --gpu 0 --log --highlyGene --epochs 500 --save-freq 1
```

## Citation
