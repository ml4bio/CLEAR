# referenced from yuqi

import scanpy as sc
import numpy as np
import pandas as pd
import os
import sys

def dropout_events(adata, drop_prob=0.0):
    adata = adata.copy()
    nnz = adata.X.nonzero()
    nnz_size = len(nnz[0])

    drop = np.random.choice(nnz_size, int(nnz_size*drop_prob))
    adata[nnz[0][drop], nnz[1][drop]] = 0

    return adata

def h5ad_to_csv(h5ad_path, save_path, filtered, normalized, drop_prob, rand_seed):
    h5ad_dir, h5ad_file = os.path.split(h5ad_path)
    h5ad_file_pre, ext = os.path.splitext(h5ad_file)
    dataset_name = h5ad_file_pre

    if os.path.exists(h5ad_path) != True:
        print("{} doesn't exist".format(h5ad_path))
        return 0
    else:
        if drop_prob == 0:
            count_csv_path = os.path.join(save_path, "{}_counts.csv".format(dataset_name))
            label_csv_path = os.path.join(save_path, "{}_labels.csv".format(dataset_name))
        else:
            count_csv_path = os.path.join(save_path, "{}_dropout{}_seed{}_counts.csv".format(dataset_name, drop_prob, rand_seed))
            label_csv_path = os.path.join(save_path, "{}_dropout{}_seed{}_labels.csv".format(dataset_name, drop_prob, rand_seed))
        if os.path.exists(count_csv_path) and os.path.exists(label_csv_path):
            print("{} already exists".format(count_csv_path))
            print("{} already exists".format(label_csv_path))
            return 0

    if os.path.exists(save_path) == False:
        os.makedirs(save_path)

    adata = sc.read_h5ad(h5ad_path)

    # filter
    if filtered == True:
        print("filtered")
        min_genes = 200  # the min_gene comes from seruat
        min_cells = 3  # the min_cell comes from seruat
        filter_mt = True
        filter_ercc = True
        # basic filtering, filter the genes and cells
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        if filter_mt:
            # filter the mitochondrial genes
            adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
            sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
            low_MT_MASK = (adata.obs.pct_counts_mt < 5)
            adata = adata[low_MT_MASK]
        if filter_ercc:
            # filter the ERCC spike-in RNAs
            adata.var['ERCC'] = adata.var_names.str.startswith('ERCC')  # annotate the group of ERCC spike-in as 'ERCC'
            sc.pp.calculate_qc_metrics(adata, qc_vars=['ERCC'], percent_top=None, log1p=False, inplace=True)
            low_ERCC_mask = (adata.obs.pct_counts_ERCC < 10)
            adata = adata[low_ERCC_mask]
    else:
        print("unfiltered")


    if drop_prob != 0:
        print("Drop Rate:{}, Seed:{}".format(drop_prob, rand_seed))
        np.random.seed(rand_seed)
        adata = dropout_events(adata, drop_prob=drop_prob)

    # save
    a = adata.to_df()
    a.to_csv(count_csv_path)
    #np.transpose(a).to_csv(count_csv_path)
    if adata.obs.get('cell_ontology_class') is not None:
        adata.obs['cell_ontology_class'].to_csv(label_csv_path)
    elif adata.obs.get('celltype') is not None:
        adata.obs['celltype'].to_csv(label_csv_path)
    elif adata.obs.get('CellType') is not None:
        adata.obs['CellType'].to_csv(label_csv_path)
    elif adata.obs.get('Group') is not None:
        adata.obs['Group'].to_csv(label_csv_path)
    else:
        adata.obs['x'].to_csv(label_csv_path)

if __name__ == '__main__':
    if len(sys.argv) == 3:
        h5ad_path = sys.argv[1]
        save_path = sys.argv[2]
        filtered = True
        normalized = False
        drop_prob = 0
        rand_seed = 0
    elif len(sys.argv) == 4:
        h5ad_path = sys.argv[1]
        save_path = sys.argv[2]
        filtered = False if sys.argv[3] != "0" else True
        normalized = 0
        drop_prob = 0
        rand_seed = 0
    elif len(sys.argv) == 5:
        h5ad_path = sys.argv[1]
        save_path = sys.argv[2]
        filtered = False if sys.argv[3] == "0" else True
        normalized = False if sys.argv[4] == "0" else True
        drop_prob = 0
        rand_seed = 0
    elif len(sys.argv) == 6:
        h5ad_path = sys.argv[1]
        save_path = sys.argv[2]
        filtered = False if sys.argv[3] == "0" else True
        normalized = False if sys.argv[4] == "0" else True
        drop_prob = float(sys.argv[5])
        rand_seed = 0
    elif len(sys.argv) == 7:
        h5ad_path = sys.argv[1]
        save_path = sys.argv[2]
        filtered = False if sys.argv[3] == "0" else True
        normalized = False if sys.argv[4] == "0" else True
        drop_prob = float(sys.argv[5])
        rand_seed = int(sys.argv[6])
    h5ad_to_csv(h5ad_path, save_path, filtered, normalized, drop_prob, rand_seed)