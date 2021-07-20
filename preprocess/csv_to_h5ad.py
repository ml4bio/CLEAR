import scanpy as sc
import os
import sys
import pandas as pd
import numpy as np

def dropout_events(adata, drop_prob=0.0):
    adata = adata.copy()
    nnz = adata.X.nonzero()
    nnz_size = len(nnz[0])

    drop = np.random.choice(nnz_size, int(nnz_size*drop_prob))
    adata[nnz[0][drop], nnz[1][drop]] = 0

    return adata

def prepocessing(csv_count_path: str,
                 save_path: str,
                 min_genes: int = 200,  # the min_gene comes from seruat
                 min_cells: int = 3,  # the min_cell comes from seruat
                 filter_mt: bool = True,
                 filter_ercc: bool = True,
                 do_cpm: bool = True,  # currently, only consider the CPM case
                 select_highly_variable: bool = False,  # optional. May be useful for some datasets
                 filtered: bool = True,
                 normalized: bool = True,
                 drop_prob: float = 0, # precentage of records that drops, choose from [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
                 rand_seed: int = 0,
                 ):
    
    # load the data into adata class
    csv_label_path = csv_count_path.replace("_counts", "_labels")
    csv_dir, csv_file = os.path.split(csv_count_path)
    csv_pre, ext = os.path.splitext(csv_file)
    dataset_name = csv_pre.replace("_counts", "")
    if os.path.exists(save_path)==False:
        os.makedirs(save_path)

    if os.path.exists(csv_count_path) != True:
        #print("Wrong file type! {}".format(csv_path))
        print("{} doesn't exist".format(csv_count_path))
        return 0
    else:
        if drop_prob == 0:
            h5ad_path = os.path.join(save_path, "{}.h5ad".format(dataset_name))
        else:
            h5ad_path = os.path.join(save_path, "{}_dropout{}_seed{}.h5ad".format(dataset_name, drop_prob, rand_seed))

        if normalized == True:
            h5ad_path = h5ad_path.replace(".h5ad", "_normalized.h5ad")

        if os.path.exists(h5ad_path):
            print("{} already exists".format(h5ad_path))
            return 0
        else:
            print("creating {}".format(h5ad_path))

    count_frame = pd.read_csv(csv_count_path, index_col=0)
    #count_frame = count_frame.reset_index(drop=True)   # reset indexes instead of the original ones
    meta_frame = pd.read_csv(csv_label_path, index_col=0)
    meta_frame.index = count_frame.index
    #meta_frame = meta_frame.reset_index(drop=True)   # reset indexes instead of the original ones
    adata = sc.AnnData(X=count_frame, obs=meta_frame)
    #adata = sc.read_csv(csv_path)
    #print(adata) #

    if filtered == True:
        print("filtered")
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

    # CJY for scVI required unnormalized data (WARNING  Make sure the registered X field in anndata contains unnormalized count data)
    if normalized == False:
        adata.write(h5ad_path)
        print(h5ad_path)
        return 0

    """
    adata_cpm = adata.copy()    # apply this to a copy so we can compare methods
    adata_cpm.raw = adata_cpm   # store a copy of the raw values before normalizing

    #sc.pp.normalize_per_cell(adata_cpm, counts_per_cell_after=1e6)   # normalize_per_cell Deprecated since version 1.3.7
    sc.pp.normalize_total(adata_cpm, target_sum=1e6)

    sc.pp.log1p(adata_cpm)
    sc.pp.scale(adata_cpm)
    h5ad_path = os.path.join(save_path, "{}_normalized.h5ad".format(dataset_name))
    adata_cpm.write(h5ad_path)
    print(h5ad_path)
        
    if select_highly_variable:
        sc.pp.highly_variable_genes(adata_cpm, min_mean=0.0125, max_mean=3, min_disp=0.5)
        adata_cpm = adata_cpm[:, adata_cpm.var.highly_variable]
        sc.pp.scale(adata, max_value=10)
        h5ad_path = os.path.join(save_path, "{}_normalized.h5ad".format(dataset_name))
        adata_cpm.write(h5ad_path)
    """

if __name__ == '__main__':
    if len(sys.argv) == 3:
        csv_file = sys.argv[1]
        save_path = sys.argv[2]
        filtered = True
        normalized = False
        drop_prob = 0
        rand_seed = 0
    elif len(sys.argv) == 4:
        csv_file = sys.argv[1]
        save_path = sys.argv[2]
        filtered = False if sys.argv[3] == "0" else True
        normalized = False
        drop_prob = 0
        rand_seed = 0
    elif len(sys.argv) == 5:
        csv_file = sys.argv[1]
        save_path = sys.argv[2]
        filtered = False if sys.argv[3] == "0" else True
        normalized = False if sys.argv[4] == "0" else True
        drop_prob = 0
        rand_seed = 0
    elif len(sys.argv) == 6:
        csv_file = sys.argv[1]
        save_path = sys.argv[2]
        filtered = False if sys.argv[3] == "0" else True
        normalized = False if sys.argv[4] == "0" else True
        drop_prob = float(sys.argv[5])
        rand_seed = 0
    elif len(sys.argv) == 7:
        csv_file = sys.argv[1]
        save_path = sys.argv[2]
        filtered = False if sys.argv[3] == "0" else True
        normalized = False if sys.argv[4] == "0" else True
        drop_prob = float(sys.argv[5])
        rand_seed = int(sys.argv[6])
    prepocessing(csv_file, save_path, filtered=filtered, normalized=normalized, drop_prob=drop_prob, rand_seed=rand_seed)
