# from h5ad to csv
target_h5ad_name_list=(
"tabula-muris-senis-facs-processed-official-annotations-Diaphragm"
"tabula-muris-senis-facs-processed-official-annotations-Limb_Muscle"
"tabula-muris-senis-facs-processed-official-annotations-Bladder"
"tabula-muris-senis-facs-processed-official-annotations-Mammary_Gland"
)
data_dir="./data/"
oh5ad_dir="${data_dir}oh5ad/"
csv_dir="${data_dir}csv/"
h5ad_dir="${data_dir}h5ad/"
h5ad_file_list=`ls $oh5ad_dir`
file_index=0
for h5ad_file in $h5ad_file_list
do
  file_index=$[ $file_index + 1 ]
  echo $file_index
  h5ad_name=$(echo ${h5ad_file/.h5ad/})
	echo $h5ad_name

  # if h5ad_name in target list
	if [ ${#target_h5ad_name_list[@]} != 0 ]   # $target_h5ad_name_list[0] != "none"
	then
	  run_flag=0
	  for target_h5ad_name in ${target_h5ad_name_list[@]}
	  do
	    if [ $target_h5ad_name == $h5ad_name ]
	    then
	      run_flag=1
	    fi
    done
    if [ $run_flag == 0 ]
    then
      echo "${h5ad_name} not in target list"
      continue
    fi
  fi

  # preprocess
	csv_path="${csv_dir}${h5ad_name}_counts.csv"
	echo $csv_path
	if [ ! -f $csv_path ]
	then
	  h5ad_path="${oh5ad_dir}${h5ad_file}"
	  python preprocess/h5ad_to_csv.py ${h5ad_path} ${csv_dir} 1 0
	  python preprocess/csv_to_h5ad.py ${csv_path} ${h5ad_dir} 0 0
	else
	  echo "${csv_path} already exists"
	fi
	#break
done