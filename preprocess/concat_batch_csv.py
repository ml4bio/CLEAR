import os
import sys
import pandas as pd

def concat_batch_csvs(root_path, dataset_name, save_path):
    dataset_path = os.path.join(root_path, dataset_name)

    csv_file_list = os.listdir(dataset_path)
    csv_file_list.sort()

    df_count_list = []
    df_label_list = []
    for csv_file in csv_file_list:
        pre, ext = os.path.splitext(csv_file)
        if ext != ".csv":
            continue
        if "celltype" in csv_file or "cell_info" in csv_file:
            df = pd.read_csv(os.path.join(dataset_path, csv_file), index_col=0)
            df_label_list.append(df)
            print("labels:{} shape:{}".format(csv_file, df.shape))
        else:
            df = pd.read_csv(os.path.join(dataset_path, csv_file), index_col=0)
            df_count_list.append(df)
            print("counts:{} shape:{}".format(csv_file, df.shape))

    save_count_path = os.path.join(save_path, "{}_counts.csv".format(dataset_name))
    df_counts = pd.concat(df_count_list, axis=0, join='inner')
    df_counts.to_csv(save_count_path)
    print("cat_counts:{} shape:{}".format(save_count_path, df_counts.shape))
    save_label_path = os.path.join(save_path, "{}_labels.csv".format(dataset_name))
    df_labels = pd.concat(df_label_list, axis=0, join='inner')
    df_labels.to_csv(save_label_path)
    print("cat_labels:{} shape:{}".format(save_label_path, df_labels.shape))

if __name__=="__main__":
    #root_path = "D:\MIP"
    #dataset_name = "dataset"
    #save_path = "D:\MIP"
    if len(sys.argv) == 4:
        root_path = sys.argv[1]
        dataset_name = sys.argv[2]
        save_path = sys.argv[3]
    if os.path.exists(save_path) == False:
        os.makedirs(save_path)

    concat_batch_csvs(root_path, dataset_name, save_path)