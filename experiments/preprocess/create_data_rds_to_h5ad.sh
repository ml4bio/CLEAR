# from rds to csv to h5ad
target_rds_name_list=()  #("baron-mouse")  #("yan" "baron-mouse" "tmuris" "lake")
rds_dir="/home/yanhan/cjy/Single-Cell-Dataset/raw_rds/"
data_dir="./data/"
ocsv_dir="${data_dir}ocsv/"
h5ad_dir="${data_dir}h5ad/"
csv_dir="${data_dir}csv/"
rds_file_list=`ls $rds_dir`
file_index=0
for rds_file in $rds_file_list
do
  file_index=$[ $file_index + 1 ]
  echo $file_index
  rds_name=$(echo ${rds_file/.rds/})
	echo $rds_name

  # if rds_name in target list
	if [ ${#target_rds_name_list[@]} != 0 ]   # $target_rds_name_list[0] != "none"
	then
	  run_flag=0
	  for target_rds_name in ${target_rds_name_list[@]}
	  do
	    if [ $target_rds_name == $rds_name ]
	    then
	      run_flag=1
	    fi
    done
    if [ $run_flag == 0 ]
    then
      echo "${rds_name} not in target list"
      continue
    fi
  fi

  # preprocess
	csv_file="${ocsv_dir}${rds_name}_counts.csv"
	echo $csv_file
	if [ ! -f $csv_file ]
	then
	  Rscript preprocess/rds_to_csv.R ${rds_dir} ${rds_name} ${ocsv_dir}
	else
	  echo "${csv_file} already exists"
	fi
	python preprocess/csv_to_h5ad.py $csv_file ${h5ad_dir} 1 0
	h5ad_file="${h5ad_dir}${rds_name}.h5ad"
  python preprocess/h5ad_to_csv.py $h5ad_file $csv_dir 0 0
	#break
done