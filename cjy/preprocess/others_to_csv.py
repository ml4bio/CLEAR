# transform dataset 1-4 batch effect to csv


import os
import sys
import pandas as pd


def txt_to_csv(txt_path, save_path, transpose=False):
    f = open(txt_path, "rb")
    a = f.readline()
    txt_dir, filename = os.path.split(txt_path)
    df = pd.read_csv(txt_path, sep="\t", index_col=0)
    save_file = os.path.join(save_path, filename.replace(".txt", ".csv"))
    if transpose == True:
        df = df.T

    df.to_csv(save_file)

    f = open(os.path.join(save_path, "csv_record.txt"), "a")
    f.write("txt_path:{}\n".format(txt_path))
    f.write("columns - len:{}, eg:{}\n".format(len(df.columns), df.columns))
    f.write("indices - len:{}, eg:{}\n".format(len(df.index), df.index))

    print("txt_path:{}".format(txt_path))
    print("columns - len:{}, eg:{}".format(len(df.columns), df.columns[0]))
    print("indices - len:{}, eg:{}".format(len(df.index), df.index[0]))

if __name__ == "__main__":
    #txt_path = r"D:\MIP\b1_exprs.txt"
    #save_path = r"D:\MIP"
    #transpose = True
    if len(sys.argv) == 4:
        txt_path = sys.argv[1]
        save_path = sys.argv[2]
        transpose = False if sys.argv[3] == "0" else True
    if os.path.exists(save_path) == False:
        os.makedirs(save_path)
    txt_to_csv(txt_path, save_path, transpose)
