import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import os

parser = argparse.ArgumentParser(description='PyTorch scRNA-seq CLR Training')

# input & ouput
parser.add_argument('--count_csv_path', type=str, default= "",
                    help='path to counts')
parser.add_argument('--label_csv_path', type=str, default= "",
                    help='path to labels')
parser.add_argument('--save_h5ad_dir', type=str, default= "",
                    help='dir to savings')
parser.add_argument('--label_colname', type=str, default="x",
                    help='column name of labels in label.csv')

# preprocessing
parser.add_argument('--CPM', action="store_true",
                    help='do count per million operation for raw counts')
parser.add_argument("--log", action="store_true",
                    help='Whether do log operation before preprocessing')
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


def preprocess_csv_to_h5ad(count_csv_path: str = "",
                           label_csv_path: str = "",
                           save_h5ad_dir: str = "",
                           label_colname: str = "x",
                           select_highly_variable_gene: bool = False,
                           do_CPM: bool = True,
                           do_log: bool = True,
                           drop_prob: float = 0.0):
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
    else:
        adata = sc.AnnData(X=count_frame)

    # adata = sc.read()

    # do preprocessing on the adata file
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

    adata = dropout_events(adata, drop_prob=drop_prob)

    if select_highly_variable_gene and not do_log:
        sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=5000)
        adata = adata[:, adata.var.highly_variable]

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

    sc.pp.scale(adata, max_value=10, zero_center=True)
    # adata[np.isnan(adata.X)] = 0
    # adata_max = np.max(adata.X)
    # adata_min = np.min(adata.X)
    # adata.X = (adata.X - adata_min)/(adata_max - adata_min)

    if save_h5ad_dir is not None:
        if os.path.exists(save_h5ad_dir) != True:
            os.makedirs(save_h5ad_dir)

        _, counts_file_name = os.path.split(count_csv_path)
        save_file_name = counts_file_name.replace(".csv", ".h5ad").replace("_counts", "")

        save_path = os.path.join(save_h5ad_dir, save_file_name)

        adata.write(save_path)
        print("successfully convert and preprocess {} to {}!".format(counts_file_name, save_file_name))

    return adata

if __name__=="__main__":
    args = parser.parse_args()

    count_csv_path = args.count_csv_path
    label_csv_path = args.label_csv_path
    save_h5ad_dir = args.save_h5ad_dir
    label_colname = args.label_colname

    processed_adata = preprocess_csv_to_h5ad(
        count_csv_path, label_csv_path, save_h5ad_dir, label_colname=label_colname,
        select_highly_variable_gene=args.highlyGene, do_CPM=args.CPM, do_log=args.log, drop_prob=args.drop_prob
    )
