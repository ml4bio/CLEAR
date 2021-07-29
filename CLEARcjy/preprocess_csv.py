import scanpy as sc


def preprocess_csv(traindir: str = "",
                       labeldir: str = "",
                       colname: str = "x",
                       select_highly_variable_gene: bool = False,
                       do_CPM: bool = True,
                       do_log: bool = True,
                       drop_prob: float = 0.0):
    # read the count matrix from the path
    count_frame = pd.read_csv(traindir, index_col=0)
    print("counts shape:{}".format(count_frame.shape))
    meta_frame = pd.read_csv(labeldir, index_col=0, header=0)
    print("labels shape:{}".format(meta_frame.shape))
    meta_frame.rename(columns={colname: 'x'}, inplace=True)
    # meta_frame.rename(columns={'cell_ontology_class': 'x'}, inplace=True)   # organ celltype
    # meta_frame.rename(columns={'CellType': 'x'}, inplace=True)   # dataset6
    # meta_frame.rename(columns={'celltype': 'x'}, inplace=True)   # dataset1
    # meta_frame.rename(columns={'Group': 'x'}, inplace=True)  # batch effect dataset3
    adata = sc.AnnData(X=count_frame, obs=meta_frame)

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

    return adata