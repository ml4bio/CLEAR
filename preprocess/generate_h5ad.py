import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import os

parser = argparse.ArgumentParser(description='PyTorch scRNA-seq CLR Training')

# input & ouput
parser.add_argument('--input_h5ad_path', type=str, default=None,
                    help='path to input h5ad files')
parser.add_argument('--input_10X_path', type=str, default=None,
                    help='path to input 10X')
parser.add_argument('--count_csv_path', type=str, default=None,
                    help='path to counts csv')
parser.add_argument('--label_csv_path', type=str, default=None,
                    help='path to labels csv')
parser.add_argument('--save_h5ad_dir', type=str, default="./",
                    help='dir to savings')

# preprocessing
parser.add_argument('--filter', action="store_false",
                    help='Whether do filtering')

parser.add_argument('--norm', action="store_false",
                    help='Whether do normalization')

parser.add_argument("--log", action="store_false",
                    help='Whether do log operation')

parser.add_argument("--scale", action="store_false",
                    help='Whether do log operation')

parser.add_argument('--CPM', action="store_true",
                    help='do count per million operation for raw counts')

parser.add_argument("--highlyGene", action="store_true",
                    help="Whether select highly variable gene")

parser.add_argument("--drop_prob", type=float, default=0.0,
                    help="simulate dropout events")


def dropout_events(adata, drop_prob=0.0):
    adata = adata.copy()
    nnz = adata.X.nonzero()
    nnz_size = len(nnz[0])

    drop = np.random.choice(nnz_size, int(nnz_size*drop_prob))
    adata[nnz[0][drop], nnz[1][drop]] = 0

    return adata


def preprocess_csv_to_h5ad(
        input_h5ad_path=None, input_10X_path=None, count_csv_path=None, label_csv_path=None, save_h5ad_dir="./",
        do_filter=False, do_norm=False, select_highly_variable_gene=False, do_log=False, do_scale=False,
        drop_prob=0.0
):
    # 1. read data from h5ad, 10X or csv files.
    if input_h5ad_path != None and input_10X_path == None and count_csv_path == None:
        adata = sc.read_h5ad(input_h5ad_path)
        print("Read data from h5ad file: {}".format(input_h5ad_path))

        _, h5ad_file_name = os.path.split(input_h5ad_path)
        save_file_name = h5ad_file_name

    elif input_10X_path != None and input_h5ad_path == None and count_csv_path == None:
        adata = sc.read_10x_mtx(input_10X_path)
        print("Read data from 10X file: {}".format(input_10X_path))

        _, input_10X_file_name = os.path.split(input_10X_path)
        save_file_name = input_10X_file_name + ".h5ad"

    elif count_csv_path != None and input_h5ad_path == None and input_10X_path == None:
        # read the count matrix from the path
        count_frame = pd.read_csv(count_csv_path, index_col=0)
        print("counts shape:{}".format(count_frame.shape))

        if label_csv_path != None:
            label_frame = pd.read_csv(label_csv_path, index_col=0, header=0)
            print("labels shape:{}".format(label_frame.shape))
            if count_frame.shape[0] != label_frame.shape[0]:
                raise Exception("The shapes of counts and labels do not match!")

            #if rename_label_colname is not None:
                #label_frame.rename(columns={label_colname: 'x'}, inplace=True)
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
        save_file_name = counts_file_name.replace(".csv", ".h5ad").replace("_counts", "")

    elif input_h5ad_path != None and count_csv_path != None:
        raise Exception("Can not address h5ad and csv files simultaneously!")


    # 2. preprocess anndata
    preprocessed_flag = do_filter | do_norm | select_highly_variable_gene | do_CPM | do_log | do_scale | drop_prob>0
    # filter operation
    if do_filter == True:
        # basic filtering, filter the genes and cells
        # sc.pp.filter_cells(adata, min_counts=5000)
        sc.pp.filter_cells(adata, min_genes=200)

        sc.pp.filter_genes(adata, min_cells=3)

        # calculate the qc metrics, then do the normalization

        # filter the mitochondrial genes
        adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        low_MT_MASK = (adata.obs.pct_counts_mt < 5)
        adata = adata[low_MT_MASK]

        # filter the ERCC spike-in RNAs
        adata.var['ERCC'] = adata.var_names.str.startswith('ERCC-')  # annotate the group of ERCC spike-in as 'ERCC'
        sc.pp.calculate_qc_metrics(adata, qc_vars=['ERCC'], percent_top=None, log1p=False, inplace=True)
        low_ERCC_mask = (adata.obs.pct_counts_ERCC < 10)
        adata = adata[low_ERCC_mask]

        preprocessed_flag = True

    # dropout operation
    if drop_prob > 0:
        adata = dropout_events(adata, drop_prob=drop_prob)

    if select_highly_variable_gene and not do_log:
        sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=5000)
        adata = adata[:, adata.var.highly_variable]

    if do_norm == True:
        sc.pp.normalize_total(adata, target_sum=1e4, exclude_highly_expressed=True)
        adata.raw = adata

    if np.max(adata.X > 100) and do_log:
        sc.pp.log1p(adata)

    # before normalization, we can select the most variant genes
    if select_highly_variable_gene and do_log:
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        adata = adata[:, adata.var.highly_variable]

    # after that, we do the linear scaling

    # sc.pp.normalize_total(adata, target_sum=1e4, exclude_highly_expressed=True)

    # log operations and scale operations will hurt the
    # contrastive between the data

    if do_scale == True:
        sc.pp.scale(adata, max_value=10, zero_center=True)
        # adata[np.isnan(adata.X)] = 0
        # adata_max = np.max(adata.X)
        # adata_min = np.min(adata.X)
        # adata.X = (adata.X - adata_min)/(adata_max - adata_min)

    # 3. save preprocessed h5ad
    if save_h5ad_dir is not None:
        if os.path.exists(save_h5ad_dir) != True:
            os.makedirs(save_h5ad_dir)

        if preprocessed_flag == True:
            save_file_name = save_file_name.replace(".h5ad", "_preprocessed.h5ad")
        save_path = os.path.join(save_h5ad_dir, save_file_name)

        adata.write(save_path)
        print("Successfully generate preprocessed file: {}".format(save_file_name))

    return adata

if __name__=="__main__":
    args = parser.parse_args()

    processed_adata = preprocess_csv_to_h5ad(
        args.input_h5ad_path, args.input_10X_path, args.count_csv_path, args.label_csv_path, args.save_h5ad_dir,
        do_filter=args.filter, do_norm=args.norm, select_highly_variable_gene=args.highlyGene,
        do_CPM=args.CPM, do_log=args.log, drop_prob=args.drop_prob, do_scale=args.scale,
    )
