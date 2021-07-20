# from h5ad to csv
target_h5ad_name_list=(
"tabula-muris-senis-facs-processed-official-annotations-Diaphragm"
"tabula-muris-senis-facs-processed-official-annotations-Limb_Muscle"
"tabula-muris-senis-facs-processed-official-annotations-Bladder"
"tabula-muris-senis-facs-processed-official-annotations-Mammary_Gland"
)
drop_prob_list=(0.1 0.3 0.6 0.8)
random_seed_list=(0 1 2 3 4)
data_dir="./data/"
# oh5ad -> csv -> h5ad
oh5ad_dir="${data_dir}oh5ad/"
csv_dir="${data_dir}csv_dropout/"
h5ad_dir="${data_dir}h5ad_dropout/"
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
  for drop_prob in ${drop_prob_list[@]}
  do
    for random_seed in ${random_seed_list[@]}
    do
      oh5ad_path="${oh5ad_dir}${h5ad_file}"
	    python preprocess/h5ad_to_csv.py ${oh5ad_path} ${csv_dir} 1 0 ${drop_prob} ${random_seed}
	    csv_path="${csv_dir}${h5ad_name}_dropout${drop_prob}_seed${random_seed}_counts.csv"
	    python preprocess/csv_to_h5ad.py ${csv_path} ${h5ad_dir} 0 0
	  done
	done
done