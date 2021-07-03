import scanpy as sc
import pandas as pd
import os
import sys

def read_csv(csv_path):
    f = open("csv_record.txt", "w")
    file_list = os.listdir(csv_path)
    file_list.sort()
    for file in file_list:
        file_path = os.path.join(csv_path, file)
        df = pd.read_csv(file_path)
        f.write("{} {} {}\n".format(file, df.shape[0], df.shape[1]))
        print("{} {} {}\n".format(file, df.shape[0], df.shape[1]))
    f.close()

def read_h5ad(h5ad_path):
    f = open("h5ad_record.txt", "w")
    file_list = os.listdir(h5ad_path)
    file_list.sort()
    for file in file_list:
        file_path = os.path.join(h5ad_path, file)
        adata = sc.read_h5ad(file_path)
        f.write("{} {} {}\n".format(file, adata.X.shape[0], adata.X.shape[1]))
        print("{} {} {}\n".format(file, adata.X.shape[0], adata.X.shape[1]))
    f.close()

if __name__ == "__main__":
    if sys.argv[1] == "csv":
        csv_path = "./data/csv_dropout"
        read_csv(csv_path)
    elif sys.argv[1] == "h5ad":
        h5ad_path = "./data/h5ad_dropout"
        read_h5ad(h5ad_path)
    else:
        raise Exception

