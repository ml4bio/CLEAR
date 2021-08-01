# encoding: utf-8
"""
@author:  Jiayang Chen
@contact: yjcmydkzgj@gmail.com

format transformation
"""

import scanpy as sc
import argparse

parser = argparse.ArgumentParser(description='PyTorch scRNA-seq CLR Training')

# input & ouput
parser.add_argument('--input_h5ad_path', type=str, default= "",
                    help='path to counts')
parser.add_argument('--count_csv_path', type=str, default= None,
                    help='path to counts')
parser.add_argument('--label_csv_path', type=str, default= None,
                    help='path to labels')
parser.add_argument('--save_h5ad_dir', type=str, default= "./",
                    help='dir to savings')
parser.add_argument('--label_colname', type=str, default="x",
                    help='column name of labels in label.csv')

def format_transformation(
        input_h5ad_path=None, input_10X_path=None, count_csv_path=None, label_csv_path=None, save_h5ad_dir="./",
):
    if input_h5ad_path != None and input_10X_path == None and count_csv_path == None:
        adata = sc.read_h5ad(input_h5ad_path)
        print("Read data from h5ad file: {}".format(input_h5ad_path))

        _, h5ad_file_name = os.path.split(input_h5ad_path)
        save_file_name = h5ad_file_name.replace(".h5ad", "_preprocessed.h5ad")

    elif input_10X_path != None and input_h5ad_path == None and count_csv_path == None:
        adata = sc.read_10x_mtx(input_10X_path)
        print("Read data from 10X file: {}".format(input_10X_path))

        _, input_10X_file_name = os.path.split(input_10X_path)
        save_file_name = input_10X_file_name + "_preprocessed.h5ad"

    elif count_csv_path != None and input_h5ad_path == None and input_10X_path == None:
        # read the count matrix from the path
        count_frame = pd.read_csv(count_csv_path, index_col=0)
        print("counts shape:{}".format(count_frame.shape))

        if label_csv_path != None:
            label_frame = pd.read_csv(label_csv_path, index_col=0, header=0)
            print("labels shape:{}".format(label_frame.shape))
            if count_frame.shape[0] != label_frame.shape[0]:
                raise Exception("The shapes of counts and labels do not match!")

            label_frame.rename(columns={label_colname: 'x'}, inplace=True)
            # label_frame.rename(columns={'cell_ontology_class': 'x'}, inplace=True)   # organ celltype
            # label_frame.rename(columns={'CellType': 'x'}, inplace=True)   # dataset6
            # label_frame.rename(columns={'celltype': 'x'}, inplace=True)   # dataset1
            # label_frame.rename(columns={'Group': 'x'}, inplace=True)  # batch effect dataset3
            label_frame.index = count_frame.index

            adata = sc.AnnData(X=count_frame, obs=label_frame)
            print("Read data from csv file: {}".format(count_csv_path))
            print("Read laebl from csv file: {}".format(label_csv_path))
        else:
            adata = sc.AnnData(X=count_frame)
            print("Read data from csv file: {}".format(count_csv_path))

        _, counts_file_name = os.path.split(count_csv_path)
        save_file_name = counts_file_name.replace(".csv", "_preprocessed.h5ad").replace("_counts", "")

    elif input_h5ad_path != None and count_csv_path != None:
        raise Exception("Can not address h5ad and csv files simultaneously!")

adata = sc.read_10x_mtx(data_path)
...(data preprocessing steps with individual paramters, see Scanpy tutorial)
adata.write('USE_FOR_CLEAR.h5ad')


#create anndata
df = pd.read_csv('count.csv', index_col=0)
cellinfo = pd.DataFrame(df.index,index=df.index,columns=['sample_index'])
geneinfo = pd.DataFrame(df.columns,index=df.columns,columns=['genes_index'])
adata = sc.AnnData(df, obs=cellinfo, var = geneinfo)

#generate h5ad file for CLEAR input
...(data preprocessing steps with individual paramters, see Scanpy tutorial)
adata.write('USE_FOR_CLEAR.h5ad')


#create anndata
df = pd.read_csv('R_FILTERED_DATA.csv', index_col=0)
cellinfo = pd.DataFrame(df.index,index=df.index,columns=['sample_index'])
geneinfo = pd.DataFrame(df.columns,index=df.columns,columns=['genes_index'])
adata = sc.AnnData(df, obs=cellinfo, var = geneinfo)

#generate h5ad file for CLEAR input
adata.write('USE_FOR_CLEAR.h5ad')