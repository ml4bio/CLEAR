
csv_dir="/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv"

raw_data_path="/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/"
raw_data6_path="${raw_data_path}dataset6"
python preprocess/txt_to_csv.py ${raw_data6_path} "${csv_dir}dataset6"



# dataset1
# /home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset1/
# /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset1
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset1/dataset1_sm_uc3.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset1" 1
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset1/sample_sm_uc3.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset1" 0

#for other methods: generate h5ad & csv
#cp /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset1/dataset1_sm_uc3.csv /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all/dataset1_counts.csv
#cp /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset1/sample_sm_uc3.csv /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all/dataset1_labels.csv
ocsv_dir="/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all/"
h5ad_dir="/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/h5ad/"
csv_dir="/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/csv/"
dataset_names=("dataset3_simul1" "dataset3_simul1_HVG" "dataset3_simul2" "dataset3_simul2_HVG" "dataset3_simul3" "dataset3_simul3_HVG" "dataset3_simul4" "dataset3_simul4_HVG" "dataset3_simul5" "dataset3_simul5_HVG" "dataset3_simul6" "dataset3_simul6_HVG")

for dataset_name in ${dataset_names[@]};
do
  ocsv_path="${ocsv_dir}${dataset_name}_counts.csv"
  h5ad_path="${h5ad_dir}${dataset_name}.h5ad"
  python preprocess/csv_to_h5ad.py ${ocsv_path} ${h5ad_dir} 1 0
  python preprocess/h5ad_to_csv.py ${h5ad_path} ${csv_dir} 0 0
done



# dataset2
# /home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset2/
# /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset2
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset2/filtered_total_batch1_seqwell_batch2_10x.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset2" 1
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset2/filtered_total_sample_ext_organ_celltype_batch.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset2" 0

# dataset3
# /home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset3/
# /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset3/simul1_dropout_005_b1_500_b2_900/counts.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul1" 0
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset3/simul1_dropout_005_b1_500_b2_900/cellinfo.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul1" 0
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset3/simul1_dropout_005_b1_500_b2_900/counts_HVG.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul1-HVG" 0
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset3/simul1_dropout_005_b1_500_b2_900/cellinfo.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul1-HVG" 0

python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset3/simul2_dropout_025_b1_500_b2_900/counts.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul2" 0
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset3/simul2_dropout_025_b1_500_b2_900/cellinfo.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul2" 0
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset3/simul2_dropout_025_b1_500_b2_900/counts_HVG.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul2-HVG" 0
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset3/simul2_dropout_025_b1_500_b2_900/cellinfo.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul2-HVG" 0

python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset3/simul3_dropout_005_b1_500_b2_450/counts.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul3" 0
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset3/simul3_dropout_005_b1_500_b2_450/cellinfo.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul3" 0
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset3/simul3_dropout_005_b1_500_b2_450/counts_HVG.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul3-HVG" 0
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset3/simul3_dropout_005_b1_500_b2_450/cellinfo.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul3-HVG" 0

python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset3/simul4_dropout_025_b1_500_b2_450/counts.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul4" 0
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset3/simul4_dropout_025_b1_500_b2_450/cellinfo.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul4" 0
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset3/simul4_dropout_025_b1_500_b2_450/counts_HVG.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul4-HVG" 0
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset3/simul4_dropout_025_b1_500_b2_450/cellinfo.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul4-HVG" 0

python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset3/simul5_dropout_005_b1_80_b2_400/counts.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul5" 0
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset3/simul5_dropout_005_b1_80_b2_400/cellinfo.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul5" 0
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset3/simul5_dropout_005_b1_80_b2_400/counts_HVG.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul5-HVG" 0
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset3/simul5_dropout_005_b1_80_b2_400/cellinfo.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul5-HVG" 0

python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset3/simul6_dropout_025_b1_80_b2_400/counts.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul6" 0
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset3/simul6_dropout_025_b1_80_b2_400/cellinfo.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul6" 0
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset3/simul6_dropout_025_b1_80_b2_400/counts_HVG.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul6-HVG" 0
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset3/simul6_dropout_025_b1_80_b2_400/cellinfo.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul6-HVG" 0

#for other methods: generate h5ad & csv
cp /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul1/counts.csv /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all/dataset3_simul1_counts.csv
cp /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul1/cellinfo.csv /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all/dataset3_simul1_labels.csv
cp /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul1-HVG/counts_HVG.csv /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all/dataset3_simul1_HVG_counts.csv
cp /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul1-HVG/cellinfo.csv /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all/dataset3_simul1_HVG_labels.csv

cp /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul2/counts.csv /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all/dataset3_simul2_counts.csv
cp /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul2/cellinfo.csv /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all/dataset3_simul2_labels.csv
cp /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul2-HVG/counts_HVG.csv /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all/dataset3_simul2_HVG_counts.csv
cp /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul2-HVG/cellinfo.csv /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all/dataset3_simul2_HVG_labels.csv

cp /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul3/counts.csv /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all/dataset3_simul3_counts.csv
cp /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul3/cellinfo.csv /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all/dataset3_simul3_labels.csv
cp /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul3-HVG/counts_HVG.csv /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all/dataset3_simul3_HVG_counts.csv
cp /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul3-HVG/cellinfo.csv /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all/dataset3_simul3_HVG_labels.csv

cp /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul4/counts.csv /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all/dataset3_simul4_counts.csv
cp /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul4/cellinfo.csv /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all/dataset3_simul4_labels.csv
cp /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul4-HVG/counts_HVG.csv /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all/dataset3_simul4_HVG_counts.csv
cp /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul4-HVG/cellinfo.csv /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all/dataset3_simul4_HVG_labels.csv

cp /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul5/counts.csv /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all/dataset3_simul5_counts.csv
cp /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul5/cellinfo.csv /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all/dataset3_simul5_labels.csv
cp /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul5-HVG/counts_HVG.csv /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all/dataset3_simul5_HVG_counts.csv
cp /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul5-HVG/cellinfo.csv /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all/dataset3_simul5_HVG_labels.csv

cp /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul6/counts.csv /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all/dataset3_simul6_counts.csv
cp /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul6/cellinfo.csv /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all/dataset3_simul6_labels.csv
cp /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul6-HVG/counts_HVG.csv /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all/dataset3_simul6_HVG_counts.csv
cp /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset3/simul6-HVG/cellinfo.csv /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all/dataset3_simul6_HVG_labels.csv

ocsv_dir="/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all/"
h5ad_dir="/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/h5ad/"
csv_dir="/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/csv/"
dataset_names=("dataset3_simul1" "dataset3_simul1_HVG" "dataset3_simul2" "dataset3_simul2_HVG" "dataset3_simul3" "dataset3_simul3_HVG" "dataset3_simul4" "dataset3_simul4_HVG" "dataset3_simul5" "dataset3_simul5_HVG" "dataset3_simul6" "dataset3_simul6_HVG")
ocsv_path="${ocsv_dir}${dataset_name}_counts.csv"
h5ad_path="${h5ad_dir}${dataset_name}.h5ad"
python preprocess/csv_to_h5ad.py ${ocsv_path} ${h5ad_dir} 1 0
python preprocess/h5ad_to_csv.py ${h5ad_path} ${csv_dir} 0 0



# dataset4
# /home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset4/
# /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset4/
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset4/myData_pancreatic_5batches.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset4" 1
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset4/mySample_pancreatic_5batches.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset4" 0


# dataset5
# /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset5
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset5/b1_exprs.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset5" 1
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset5/b1_celltype.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset5" 0
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset5/b2_exprs.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset5" 1
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset5/b2_celltype.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset5" 0

python preprocess/concat_batch_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv" "dataset5" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all"


# dataset6
# /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset6
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset6/b1_exprs.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset6" 1
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset6/b1_celltype.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset6" 0
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset6/b2_exprs.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset6" 1
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset6/b2_celltype.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset6" 0
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset6/b3_exprs.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset6" 1
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset6/b3_celltype.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset6" 0

python preprocess/concat_batch_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv" "dataset6" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all"

# dataset7
# /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset7
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset7/b1_exprs.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset7" 1
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset7/b1_celltype.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset7" 0
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset7/b2_exprs.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset7" 1
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset7/b2_celltype.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset7" 0

python preprocess/concat_batch_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv" "dataset7" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all"

# dataset8
# /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset8


# dataset10
# /home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset7
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset10/b1_exprs.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset10" 1
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset10/b1_celltype.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset10" 0
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset10/b2_exprs.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset10" 1
python preprocess/txt_to_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/batch-effect-removal/batch_effect/dataset10/b2_celltype.txt" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/dataset10" 0

python preprocess/concat_batch_csv.py "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv" "dataset10" "/home/yanhan/cjy/Single-Cell-Dataset/Scarlet/data/ocsv/all"
