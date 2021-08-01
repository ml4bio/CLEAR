import ItClust as ic
import scanpy.api as sc
import os
from numpy.random import seed
from tensorflow import set_random_seed
import numpy as np
import warnings
import sys

import metrics as M

os.environ["CUDA_VISIBLE_DEVICES"]="1"
warnings.filterwarnings("ignore")
#import sys
#!{sys.executable} -m pip install 'scanpy==1.4.4.post1'
#Set seeds
seed(20180806)
np.random.seed(10)
set_random_seed(20180806)   # on GPU may be some other default

def source_index_map(x):
    return "source_"+x

def target_index_map(x):
    return x.replace("_", "").replace("-", "")#"target_"+x

def train(data_path, dataset_name, save_path):
    """
    :param data_path: e.g. data/    including csv h5ad
    :param dataset_name: "baron-human"
    :param save_path: e.g. results/
    """
    # load data
    #adata_source = sc.read("./data/pancreas/Bh.h5ad")   # must keep index of source and target different
    #adata_target = sc.read("./data/pancreas/smartseq2.h5ad")

    #"""
    source_dataset_name = "baron-human"
    source_path_to_anndata = os.path.join(data_path, "{}.h5ad".format(source_dataset_name))
    if os.path.exists(source_path_to_anndata)==False:
        raise Exception("{} doesn't exist".format(source_path_to_anndata))
    adata_source = sc.read_h5ad(source_path_to_anndata)
    adata_source.obs["celltype"] = adata_source.obs["x"] if adata_source.obs.get("x") is not None else adata_source.obs["cell_ontology_class"]

    #new_obs = adata_source.obs.loc[:, ["celltype"]]
    #adata_source.obs = new_obs
    #adata_source.var = pd.DataFrame(index=adata_source.var.index)
    # adata_source.obs.rename(index=source_index_map, inplace=True)   # Invalid
    #adata_source.obs = adata_source.obs.reset_index(drop=True)

    #adata_source.var = adata_source.var.reset_index(drop=True)

    #adata_source.var.index = adata_source.var.index.tolist()
    #"""

    #"""
    target_path_to_anndata = os.path.join(data_path, "{}.h5ad".format(dataset_name))
    if os.path.exists(target_path_to_anndata)==False:
        raise Exception("{} doesn't exist".format(target_path_to_anndata))
    adata_target = sc.read_h5ad(target_path_to_anndata)
    if adata_target.obs.get("x") is not None:
        adata_target.obs["celltype"] = adata_target.obs["x"]
    elif adata_target.obs.get("cell_ontology_class") is not None:
        adata_target.obs["celltype"] = adata_target.obs["cell_ontology_class"]
    elif adata_target.obs.get("CellType") is not None:
        adata_target.obs["celltype"] = adata_target.obs["CellType"]
    elif adata_target.obs.get("Group") is not None:
        adata_target.obs["celltype"] = adata_target.obs["Group"]
    #adata_target.obs.index = adata_target.obs.index.tolist()
    #new_obs = adata_target.obs.loc[:, ["celltype"]]
    #adata_target.obs = new_obs
    #adata_target.var = pd.DataFrame(index=adata_target.var.index)
    #adata_target.var.rename(index=target_index_map, inplace=True)  # Invalid
    #adata_target.obs = adata_target.obs.reset_index(drop=True)

    #adata_target.var = adata_target.var.reset_index(drop=True)

    #adata_target.var.index = adata_target.var.index.tolist()
    #"""

    # train model
    clf = ic.transfer_learning_clf()
    clf.fit(adata_source, adata_target)

    # prediction
    pred, prob, celltype_pred = clf.predict(write=False)
    #pred.head()

    # compute metrics
    labels = adata_target.obs["celltype"]
    labels_pred = pred.cluster.values
    #ari_score = ARI(labels, labels_pred)
    #nmi_score = NMI(labels, labels_pred)
    #print("\nARI: %.4f\nNMI: %.4f"  % (ari_score, nmi_score))
    metrics = M.compute_metrics(labels, labels_pred)
    print("{}: \n {}".format(dataset_name, metrics))

    # save txt
    save_path = os.path.join(save_path, "ItClust")
    if os.path.exists(save_path) != True:
        os.makedirs(save_path)
    txt_path = os.path.join(save_path, "metric_ItClust.txt")
    f = open(txt_path, "a")
    #f.write("{} {} {}\n".format(dataset_name, ari_score, nmi_score))  #:.4f
    record_string = dataset_name
    for key in metrics.keys():
        record_string += " {}".format(metrics[key])
    record_string += "\n"
    f.write(record_string)
    f.close()

    # save data
    #adata.write_h5ad(os.path.join(save_path, "scVI_{}.h5ad".format(dataset_name)))
    np.savetxt(os.path.join(save_path, "feature_ItClust_{}.csv".format(dataset_name)), clf.adata_test.obsm["X_Embeded_z" + str(clf.save_atr)], delimiter=',')
    adata_target.obs["celltype"].to_csv(os.path.join(save_path, "gt_label_ItClust_{}.csv".format(dataset_name)))
    pred.cluster.to_csv(os.path.join(save_path, "pd_label_ItClust_{}.csv".format(dataset_name)))
    #np.save(os.path.join(save_path, "ItClust_{}_X_Embeded.npy".format(dataset_name)), clf.adata_test.obsm["X_Embeded_z" + str(clf.save_atr)])
    print("Successfully save data")

if __name__ == "__main__":
    if len(sys.argv) == 4:
        root_path = sys.argv[1]
        dataset_name = sys.argv[2]
        save_path = sys.argv[3]
        train(root_path, dataset_name, save_path)
    else:
        raise Exception("Wrong Argv Num")


