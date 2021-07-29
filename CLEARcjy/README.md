This is a PyTorch implementation of MoCo for single cell data analysis


To run the code for unsupervised learning:

python main_pcl.py --lr 1 --batch-size 512 --pcl-r 1024 --cos --count_data "/path/to/count/" --label_data "/path/to/label" --log --highlyGene --exp-dir exp


the framework of CLEAR includes three parts.

1.preprocess the input dataset

2.training the encoder

3.compute metrics
