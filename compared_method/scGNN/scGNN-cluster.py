# https://github.com/juexinwang/scGNN.git
import pandas as pd
from sklearn.metrics import adjusted_rand_score as ARI
from sklearn.metrics import normalized_mutual_info_score as NMI
import os
import sys

def train(root_path, dataset_name, save_path):
    labels_path = os.path.join(root_path, dataset_name + "_labels.csv")

    save_path = os.path.join(save_path, "scGNN/")

    parent_path = os.path.split(root_path.strip("/"))[0]
    print(parent_path)
    if "dropout" in dataset_name:
        ltmg_path = os.path.join(parent_path, "ltmg_dropout")
    else:
        ltmg_path = os.path.join(parent_path, "ltmg")
    ltmg_path = ltmg_path + "/"
    print(ltmg_path)

    if os.path.exists(os.path.join(ltmg_path, dataset_name)) != True:
        cmd_line1 = "python -W ignore method/scGNN/PreprocessingscGNN.py --datasetName {} --datasetDir {}" \
                    " --LTMGDir {} --filetype CSV --geneSelectnum 2000 --inferLTMGTag --transpose".format(dataset_name, root_path, ltmg_path)
        os.system(cmd_line1)
    else:
        print("Already Preprocessed")

    cmd_line2 = "python -W ignore method/scGNN/scGNN.py --datasetName {} --datasetDir {} --LTMGDir {} " \
                "--outputDir {} --EM-iteration 2 --Regu-epochs 50 --EM-epochs 20 --quickmode --nonsparseMode " \
                "--regulized-type LTMG".format(dataset_name, ltmg_path, ltmg_path, save_path)
    os.system(cmd_line2)

    gt_df = pd.read_csv(labels_path)
    label_gt = gt_df["x"].values if gt_df.get("x") is not None else gt_df["cell_ontology_class"].values


    if os.path.exists(save_path) != True:
        os.makedirs(save_path)

    pd_path = os.path.join(save_path, "{}_results.txt".format(dataset_name))
    pd_df = pd.read_csv(pd_path)
    label_pd = pd_df['Celltype']

    print(label_gt.shape)
    print(label_pd.shape)

    ari_score = ARI(label_gt, label_pd)
    nmi_score = NMI(label_gt, label_pd)
    print("\nARI: %.4f\nNMI: %.4f" % (ari_score, nmi_score))

    # save txt
    txt_path = os.path.join(save_path, "metric_scGNN.txt")
    f = open(txt_path, "a")
    f.write("{} {} {}\n".format(dataset_name, ari_score, nmi_score))  #:.4f
    f.close()

if __name__ == "__main__":
    if len(sys.argv) == 4:
        root_path = sys.argv[1]
        dataset_name = sys.argv[2]
        save_path = sys.argv[3]
        train(root_path, dataset_name, save_path)
    else:
        raise Exception("Wrong Argv Num")