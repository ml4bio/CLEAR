from sklearn.metrics import homogeneity_score
from sklearn.metrics import completeness_score
from sklearn.metrics import v_measure_score
from sklearn.metrics import adjusted_rand_score as ARI
from sklearn.metrics import normalized_mutual_info_score as NMI

import numpy as np
import pandas as pd
import os
import sys
import csv

def traverse_folder_compute_metrics(root_path, save_path):
    for filename in os.listdir(root_path):
        if "gt_label" not in filename:
            continue
        pre_filename, ext = os.path.splitext(filename)
        _, _, method_name, dataset_name = pre_filename.split("_", 3)

        gt_label_path = os.path.join(root_path, filename)
        pd_label_path = os.path.join(root_path, filename.replace("gt_label", "pd_label"))
        gt_label = pd.read_csv(gt_label_path, index_col=0).iloc[:, 0]
        pd_label = pd.read_csv(pd_label_path, index_col=0).iloc[:, 0]
        print(gt_label.shape)
        print(pd_label.shape)

        metrics = compute_metrics(gt_label, pd_label)

        print("{} - {}: \n {}".format(dataset_name, method_name, metrics))

        # save txt
        #save_path = os.path.join(save_path, method_name)
        if os.path.exists(save_path) != True:
            os.makedirs(save_path)
        #txt_path = os.path.join(save_path, "new_metrics.txt")
        #f = open(txt_path, "a")
        # f.write("{} {} {}\n".format(dataset_name, ari_score, nmi_score))

        #record_string = dataset_name + " " + method_name
        #for key in metrics.keys():
        #    record_string += " {}".format(metrics[key])
        #record_string += "\n"
        #f.write(record_string)
        #f.close()

        csv_path = os.path.join(save_path, "new_metrics.csv")
        need_head = False
        if os.path.exists(csv_path) != True:
            need_head = True
        f = open(csv_path, "a", encoding="utf-8’", newline="")
        csv_writer = csv.writer(f)
        record = [dataset_name, method_name]
        head = ["dataset", "method"]
        for key in metrics.keys():
            record.append(metrics[key])
            head.append(key)
        if need_head == True:
            csv_writer.writerow(head)
        csv_writer.writerow(record)
        f.close()



def compute_metrics(y_true, y_pred):
    metrics = {}
    metrics["ARI"] = ARI(y_true, y_pred)
    metrics["NMI"] = NMI(y_true, y_pred)
    metrics["CA"] = cluster_acc(y_true, y_pred)
    metrics["JI"] = Jaccard_index(y_true, y_pred)
    metrics["CS"] = completeness_score(y_true, y_pred)
    metrics["HS"] = homogeneity_score(y_true, y_pred)
    metrics["VMS"] = v_measure_score(y_true, y_pred)

    return metrics


# additional metrics for clustering.
# 1. Cluster accuracy (CA)
def cluster_acc(y_true, y_pred):
    """
    Calculate clustering accuracy. Require scikit-learn installed
    # Arguments
        y: true labels, numpy.array with shape `(n_samples,)`
        y_pred: predicted labels, numpy.array with shape `(n_samples,)`
    # Return
        accuracy, in [0,1]
    """
    if isinstance(y_true, pd.DataFrame):
        y_true = y_true.values

    if isinstance(y_pred, pd.DataFrame):
        y_pred = y_pred.values

    # CJY 2021.6.20 for batch effect dataset
    #y_true = y_true.to_list()
    y_pred = y_pred.astype(np.int64)

    label_to_number = {label: number for number, label in enumerate(set(y_true))}
    label_numerical = np.array([label_to_number[i] for i in y_true])

    y_true = label_numerical.astype(np.int64)
    assert y_pred.size == y_true.size
    D = max(y_pred.max(), y_true.max()) + 1
    w = np.zeros((D, D), dtype=np.int64)
    for i in range(y_pred.size):
        w[y_pred[i], y_true[i]] += 1
    # from sklearn.utils.linear_assignment_ import linear_assignment
    # ind = linear_assignment(w.max() - w)
    # return sum([w[i, j] for i, j in ind]) * 1.0 / y_pred.size
    # https://stackoverflow.com/questions/62390517/no-module-named-sklearn-utils-linear-assignment
    from scipy.optimize import linear_sum_assignment as linear_assignment
    row_ind, col_ind = linear_assignment(w.max() - w)
    return sum([w[i, j] for i, j in zip(row_ind, col_ind)]) * 1.0 / y_pred.size


# 2. Jaccard score (JS)
# Please refer to https://scikit-learn.org/stable/modules/generated/sklearn.metrics.jaccard_score.html
def Jaccard_index(y_true, y_pred):
    from sklearn.metrics.cluster import pair_confusion_matrix
    contingency = pair_confusion_matrix(y_true, y_pred)
    JI = contingency[1,1]/(contingency[1,1]+contingency[0,1]+contingency[1,0])
    return JI

# 3. Homogeneity score.
# https://scikit-learn.org/stable/modules/generated/sklearn.metrics.homogeneity_score.html#sklearn.metrics.homogeneity_score
# Note that Homogeneity score is asysmmetric

# 4. Completeness score 
# https://scikit-learn.org/stable/modules/generated/sklearn.metrics.completeness_score.html#sklearn.metrics.completeness_score
# It's also asysmmetric

# 5. V measure score
# https://scikit-learn.org/stable/modules/generated/sklearn.metrics.v_measure_score.html#sklearn.metrics.v_measure_score


if __name__=="__main__":
    #root_path = "/home/yanhan/cjy/Single-Cell-Dataset/Single-Cell-Cluster/result/scVI"
    #save_path = "/home/yanhan/cjy/Single-Cell-Dataset/Single-Cell-Cluster/result/"
    #root_path = r"C:\Users\cjy\Desktop\Seurat"
    #save_path = r"C:\Users\cjy\Desktop\"
    if len(sys.argv) == 3:
        root_path = sys.argv[1]
        save_path = sys.argv[2]
        traverse_folder_compute_metrics(root_path, save_path)
    else:
        raise Exception("Wrong Argv Num")