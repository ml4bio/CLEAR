import os
import sys
import scanpy as sc
import numpy as np
from sklearn.cluster import KMeans

import metrics as M


#from .metrics import compute_metrics


def train(data_path, dataset_name, save_path):
    """
    :param data_path: e.g. data/    including csv h5ad
    :param dataset_name: "baron-human"
    :param save_path: e.g. results/
    """
    # load data
    path_to_anndata = os.path.join(data_path, "{}.h5ad".format(dataset_name))
    if os.path.exists(path_to_anndata)==False:
        raise Exception("{} doesn't exist".format(path_to_anndata))
    adata = sc.read_h5ad(path_to_anndata)
    #path_to_label = os.path.join(data_path, "csv", "{}_labels.csv".format(dataset_name))
    #if os.path.exists(path_to_label)==False:
    #    raise Exception("{} doesn't exist".format(path_to_label))
    #alabel = pd.read_csv(path_to_label)
    print(adata)
    if adata.obs.get("cell_ontology_class") is not None:
        labels = adata.obs["cell_ontology_class"]
    elif adata.obs.get("x") is not None:
        labels = adata.obs["x"]
    elif adata.obs.get("celltype") is not None:
        labels = adata.obs["celltype"]
    elif adata.obs.get("CellType") is not None:
        labels = adata.obs["CellType"]
    elif adata.obs.get("Group") is not None:
        labels = adata.obs["Group"]

    # cluster
    k = labels.unique().shape[0]
    print("num cluster class: {}".format(k))
    kmeans = KMeans(n_clusters=k, random_state=0).fit(adata.X)
    adata.obs['kmeans'] = kmeans.labels_.astype(str)

    # compute metrics
    labels_pred = adata.obs['kmeans']  # ["leiden_totalVI"]

    #ari_score = ARI(labels, labels_pred)
    #nmi_score = NMI(labels, labels_pred)
    #print("\nARI: %.4f\nNMI: %.4f" % (ari_score, nmi_score))
    metrics = M.compute_metrics(labels, labels_pred)
    print("{}: \n {}".format(dataset_name, metrics))


    # save txt
    save_path = os.path.join(save_path, "DirectKMeans")
    if os.path.exists(save_path)!=True:
        os.makedirs(save_path)
    txt_path = os.path.join(save_path, "metric_DirectKMeans.txt")
    f = open(txt_path, "a")
    #f.write("{} {} {}\n".format(dataset_name, ari_score, nmi_score))
    record_string = dataset_name
    for key in metrics.keys():
        record_string += " {}".format(metrics[key])
    record_string += "\n"
    f.write(record_string)
    f.close()

    # save data
    #adata.write_h5ad(os.path.join(save_path, "scVI_{}.h5ad".format(dataset_name)))
    np.savetxt(os.path.join(save_path, "feature_DirectKMeans_{}.csv".format(dataset_name)), adata.X, delimiter=',')
    labels.to_csv(os.path.join(save_path, "gt_label_DirectKMeans_{}.csv".format(dataset_name)))
    labels_pred.to_csv(os.path.join(save_path, "pd_label_DirectKMeans_{}.csv".format(dataset_name)))
    #np.save(os.path.join(save_path, "scVI_{}_X_scVI.npy".format(dataset_name)), adata.obsm["X_scVI"])
    #np.save(os.path.join(save_path, "scVI_{}_X_normalized_scVI.npy".format(dataset_name)), adata.obsm["X_normalized_scVI"].values)
    #adata.obsm["X_normalized_scVI"].to_csv(os.path.join(save_path, "scVI_{}_X_normalized_scVI.csv".format(dataset_name)))
    print("Successfully save data")


if __name__ == "__main__":
    if len(sys.argv) == 4:
        root_path = sys.argv[1]
        dataset_name = sys.argv[2]
        save_path = sys.argv[3]
        train(root_path, dataset_name, save_path)
    else:
        raise Exception("Wrong Argv Num")