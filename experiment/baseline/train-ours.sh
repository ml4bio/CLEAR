target_h5ad_name_list=(
#"baron-mouse" "deng" "hrvatin" "kolodziejczyk" "muraro" "pollen" "yan"
#"tabula-muris-senis-facs-processed-official-annotations-Diaphragm"
"tabula-muris-senis-facs-processed-official-annotations-Limb_Muscle"
#"tabula-muris-senis-facs-processed-official-annotations-Bladder"
#"tabula-muris-senis-facs-processed-official-annotations-Mammary_Gland"
)  #("yan" "baron-mouse" "tmuris" "lake")  #"campbell"
#rds_dir="/home/yanhan/cjy/Single-Cell-Dataset/raw_rds/"
data_dir="./data/"
csv_dir="${data_dir}csv/"
h5ad_dir="${data_dir}h5ad/"
save_dir="./result/"
h5ad_file_list=`ls $h5ad_dir`
file_index=0
for h5ad_file in $h5ad_file_list
do
  file_index=$[ $file_index + 1 ]
  echo $file_index
  dataset_name=$(echo ${h5ad_file/.h5ad/})
	echo $dataset_name

	# if rds_name in target list
	if [ ${#target_h5ad_name_list[@]} != 0 ]
	then
	  run_flag=0
	  for target_h5ad_name in ${target_h5ad_name_list[@]}
	  do
	    if [ $target_h5ad_name == $dataset_name ]
	    then
	      run_flag=1
	    fi
    done
    if [ $run_flag == 0 ]
    then
      echo "${dataset_name} not in target list"
      continue
    fi
  fi

	#python method/scVI-train.py ${h5ad_dir} ${dataset_name} ${save_dir}
	#Rscript method/scDHA-train.R ${csv_dir} ${dataset_name} ${save_dir}
	#python method/ItClust-train.py ${h5ad_dir} ${dataset_name} ${save_dir}

	if [ $dataset_name == "baron-mouse" ]
	then
	  python method/Scarlet/main_pcl_cjy.py --lr 10 --batch-size 512 --pcl-r 1024 --cos --count_data "${csv_dir}${dataset_name}_counts.csv" --label_data "${csv_dir}${dataset_name}_labels.csv" --exp-dir exp_tmp --gpu 0 --epochs 500
	elif [ $dataset_name == "deng" ]
	then
	  python method/Scarlet/main_pcl_cjy.py --lr 1 --batch-size 512 --pcl-r 1024 --cos --count_data "${csv_dir}${dataset_name}_counts.csv" --label_data "${csv_dir}${dataset_name}_labels.csv" --exp-dir exp_tmp --gpu 0 --log --epochs 500
	elif [ $dataset_name == "hrvatin" ]
	then
	  python method/Scarlet/main_pcl_cjy.py --lr 1 --batch-size 512 --pcl-r 1024 --cos --count_data "${csv_dir}${dataset_name}_counts.csv" --label_data "${csv_dir}${dataset_name}_labels.csv" --exp-dir exp_tmp --gpu 0 --log --highlyGene --epochs 500
	elif [ $dataset_name == "kolodziejczyk" ]
	then
	  python method/Scarlet/main_pcl_cjy.py --lr 1 --batch-size 512 --pcl-r 1024 --cos --count_data "${csv_dir}${dataset_name}_counts.csv" --label_data "${csv_dir}${dataset_name}_labels.csv" --exp-dir exp_tmp --gpu 0 --highlyGene --epochs 500
	elif [ $dataset_name == "muraro" ]
	then
	  python method/Scarlet/main_pcl_cjy.py --lr 5e-1 --batch-size 512 --pcl-r 1024 --cos --count_data "${csv_dir}${dataset_name}_counts.csv" --label_data "${csv_dir}${dataset_name}_labels.csv" --exp-dir exp_tmp --gpu 0 --log --highlyGene --epochs 500
	elif [ $dataset_name == "pollen" ]
	then
	  python method/Scarlet/main_pcl_cjy.py --lr 1 --batch-size 512 --pcl-r 1024 --cos --count_data "${csv_dir}${dataset_name}_counts.csv" --label_data "${csv_dir}${dataset_name}_labels.csv" --exp-dir exp_tmp --gpu 0 --log --highlyGene --epochs 500
	elif [ $dataset_name == "yan" ]
	then
	  python method/Scarlet/main_pcl_cjy.py --lr 1 --batch-size 512 --pcl-r 1024 --cos --count_data "${csv_dir}${dataset_name}_counts.csv" --label_data "${csv_dir}${dataset_name}_labels.csv" --exp-dir exp_tmp --gpu 0 --log --highlyGene --epochs 500
	elif [ $dataset_name == "tabula-muris-senis-facs-processed-official-annotations-Diaphragm" ]
	then
	  python method/Scarlet/main_pcl_cjy.py --lr 5e-1 --batch-size 512 --pcl-r 1024 --cos --count_data "${csv_dir}${dataset_name}_counts.csv" --label_data "${csv_dir}${dataset_name}_labels.csv" --exp-dir exp_tmp --gpu 0 --log --highlyGene --epochs 500 --save-freq 1
	elif [ $dataset_name == "tabula-muris-senis-facs-processed-official-annotations-Limb_Muscle" ]
	then
	  python method/Scarlet/main_pcl_cjy.py --lr 1 --batch-size 512 --pcl-r 1024 --cos --count_data "${csv_dir}${dataset_name}_counts.csv" --label_data "${csv_dir}${dataset_name}_labels.csv" --exp-dir exp_tmp --gpu 0 --epochs 500 --save-freq 1
	elif [ $dataset_name == "tabula-muris-senis-facs-processed-official-annotations-Bladder" ]
	then
	  python method/Scarlet/main_pcl_cjy.py --lr 5 --batch-size 512 --pcl-r 1024 --cos --count_data "${csv_dir}${dataset_name}_counts.csv" --label_data "${csv_dir}${dataset_name}_labels.csv" --exp-dir exp_tmp --gpu 0 --log --highlyGene --epochs 500 --save-freq 1
	elif [ $dataset_name == "tabula-muris-senis-facs-processed-official-annotations-Mammary_Gland" ]
	then
	  python method/Scarlet/main_pcl_cjy.py --lr 5e-1 --batch-size 512 --pcl-r 1024 --cos --count_data "${csv_dir}${dataset_name}_counts.csv" --label_data "${csv_dir}${dataset_name}_labels.csv" --exp-dir exp_tmp --gpu 0 --log --highlyGene --epochs 500 --save-freq 1
  else
   echo "No Match DataName"
  fi
	#break
done