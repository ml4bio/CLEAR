# only use dataset1 and dataset6
ocsv_dir="/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all/"
h5ad_dir="/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/h5ad/"
csv_dir="/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/csv/"
ocsv_path="${ocsv_dir}dataset1_counts.csv"
h5ad_path="${h5ad_dir}dataset1.h5ad"
python preprocess/csv_to_h5ad.py ${ocsv_path} ${h5ad_dir} 1 0
python preprocess/h5ad_to_csv.py ${h5ad_path} ${csv_dir} 0 0

ocsv_path="${ocsv_dir}dataset6_counts.csv"
h5ad_path="${h5ad_dir}dataset6.h5ad"
python preprocess/csv_to_h5ad.py ${ocsv_path} ${h5ad_dir} 1 0
python preprocess/h5ad_to_csv.py ${h5ad_path} ${csv_dir} 0 0