from PIL import ImageFilter
import random
from anndata._core.anndata import AnnData
import torchvision.datasets as datasets
from torch.utils.data import Dataset
import scanpy as sc
import torch
from copy import deepcopy
import numpy as np
import pandas as pd



class TwoCropsTransform:
    """Take two random crops of one image as the query and key."""

    def __init__(self, base_transform):
        self.base_transform = base_transform

    def __call__(self, x):
        q = self.base_transform(x)
        k = self.base_transform(x)
        return [q, k]


class GaussianBlur(object):
    """Gaussian blur augmentation in SimCLR https://arxiv.org/abs/2002.05709"""

    def __init__(self, sigma=[.1, 2.]):
        self.sigma = sigma

    def __call__(self, x):
        sigma = random.uniform(self.sigma[0], self.sigma[1])
        x = x.filter(ImageFilter.GaussianBlur(radius=sigma))
        return x
    

class ImageFolderInstance(datasets.ImageFolder):
    def __getitem__(self, index):
        path, target = self.samples[index]
        sample = self.loader(path)
        if self.transform is not None:
            sample = self.transform(sample)           
        return sample, index





class scRNAMatrixInstance(Dataset):
    def __init__(self,
                 adata: AnnData = None,
                 adata_sct: AnnData = None,
                 obs_label_colname: str = "x",
                 transform: bool = False,
                 args_transformation: dict = {}
                 ):

        super().__init__()

        self.adata = adata

        # data
        # scipy.sparse.csr.csr_matrix or numpy.ndarray
        if isinstance(self.adata.X, np.ndarray):
            self.data = self.adata.X
        else:
            self.data = self.adata.X.toarray()

        # label (if exist, build the label encoder)
        if self.adata.obs.get(obs_label_colname) is not None:
            self.label = self.adata.obs[obs_label_colname]
            self.unique_label = list(set(self.label))
            self.label_encoder = {k: v for k, v in zip(self.unique_label, range(len(self.unique_label)))}
            self.label_decoder = {v: k for k, v in self.label_encoder.items()}
        else:
            self.label = None
            print("Can not find corresponding labels")

        # do the transformation
        self.transform = transform
        self.num_cells, self.num_genes = self.adata.shape
        self.args_transformation = args_transformation
        
        self.dataset_for_transform = deepcopy(self.data)
        #import pdb; pdb.set_trace()
        if isinstance(adata_sct.X, np.ndarray):
            adata_sct = adata_sct.X
        else:
            adata_sct = adata_sct.X.toarray()
        self.sct_dataset_for_transform  = deepcopy(adata_sct)
        
    def RandomTransform(self, sample, index):
        #tr_sample = deepcopy(sample)
        tr = transformation(self.dataset_for_transform,self.sct_dataset_for_transform[index], sample)
        
        # the crop operation

        # Mask
        tr.random_mask(self.args_transformation['mask_percentage'], self.args_transformation['apply_mask_prob'])

        # (Add) Gaussian noise
        tr.random_gaussian_noise(self.args_transformation['noise_percentage'], self.args_transformation['sigma'], self.args_transformation['apply_noise_prob'])
        # swap with NB regressed mu
        tr.random_nb_swap(self.args_transformation['nb_percentage'], self.args_transformation['apply_nb_prob'])
        # inner swap
        tr.random_swap(self.args_transformation['swap_percentage'], self.args_transformation['apply_swap_prob'])
        
        # cross over with one instance
        tr.instance_crossover(self.args_transformation['cross_percentage'], self.args_transformation['apply_cross_prob'])

        # cross over with many instances
        tr.tf_idf_based_replacement(self.args_transformation['change_percentage'], self.args_transformation['apply_mutation_prob'], True)
        tr.ToTensor()

        return tr.cell_profile


    def __getitem__(self, index):
        
        sample = self.data[index]

        if self.label is not None:
            label = self.label_encoder[self.label[index]]
        else:
            label = -1

        if self.transform:
            sample_1 = self.RandomTransform(sample, index)
            sample_2 = self.RandomTransform(sample, index)
            sample = [sample_1, sample_2]
        
        return sample, index, label

    def __len__(self):
        return self.adata.X.shape[0]


class transformation():
    
    def __init__(self, 
                 dataset,
                 sct_profile,
                 cell_profile):
        self.dataset = dataset
        self.sct_profile = sct_profile
        self.cell_profile = deepcopy(cell_profile)
        self.gene_num = len(self.cell_profile)
        self.cell_num = len(self.dataset)
    
    
    def build_mask(self, masked_percentage: float):
        mask = np.concatenate([np.ones(int(self.gene_num * masked_percentage), dtype=bool), 
                               np.zeros(self.gene_num - int(self.gene_num * masked_percentage), dtype=bool)])
        np.random.shuffle(mask)
        return mask
    

    def RandomCrop(self,
                   crop_percentage=0.8):
        mask = self.build_mask(crop_percentage)
        self.cell_profile = self.cell_profile[mask]
        self.gene_num = len(self.cell_profile)
        self.dataset = self.dataset[:,mask]



    def random_mask(self, 
                    mask_percentage: float = 0.15, 
                    apply_mask_prob: float = 0.5):

        s = np.random.uniform(0,1)
        if s<apply_mask_prob:
            # create the mask for mutation
            mask = self.build_mask(mask_percentage)
            
            # do the mutation with prob
        
            self.cell_profile[mask] = 0



    def random_gaussian_noise(self, 
                              noise_percentage: float=0.2, 
                              sigma: float=0.5, 
                              apply_noise_prob: float=0.3):

        s = np.random.uniform(0,1)
        if s<apply_noise_prob:
            # create the mask for mutation
            mask = self.build_mask(noise_percentage)
            
            # create the noise
            noise = np.random.normal(0, sigma, int(self.gene_num*noise_percentage))
            
            # do the mutation (maybe not add, simply change the value?)
            self.cell_profile[mask] += noise
    
    def random_nb_swap(self,
                       nb_swap_percentage: float=0.8,
                       apply_nb_prob: float=0.5):
            s = np.random.uniform(0,1)
            if s<apply_nb_prob:
                mask = self.build_mask(nb_swap_percentage)
                self.sct_profile[mask], self.cell_profile[mask]  = self.cell_profile[mask], self.sct_profile[mask]

    def random_swap(self,
                    swap_percentage: float=0.1,
                    apply_swap_prob: float=0.5):

        ##### for debug
        #     from copy import deepcopy
        #     before_swap = deepcopy(cell_profile)
        s = np.random.uniform(0,1)
        if s<apply_swap_prob:
            # create the number of pairs for swapping 
            swap_instances = int(self.gene_num*swap_percentage/2)
            swap_pair = np.random.randint(self.gene_num, size=(swap_instances,2))
            
            # do the inner crossover with p
        
            self.cell_profile[swap_pair[:,0]], self.cell_profile[swap_pair[:,1]] = \
                self.cell_profile[swap_pair[:,1]], self.cell_profile[swap_pair[:, 0]]



    def instance_crossover(self,
                           cross_percentage: float=0.25,
                           apply_cross_prob: float=0.4):
        
        # it's better to choose a similar profile to crossover
        
        s = np.random.uniform(0,1)
        if s<apply_cross_prob:
            # choose one instance for crossover
            cross_idx = np.random.randint(self.cell_num)
            cross_instance = self.dataset[cross_idx]
            
            # build the mask
            mask = self.build_mask(cross_percentage)
            
            # apply instance crossover with p
            
            tmp = cross_instance[mask].copy()
        
            cross_instance[mask], self.cell_profile[mask]  = self.cell_profile[mask], tmp


    def tf_idf_based_replacement(self, 
                                 change_percentage: float=0.25,
                                 apply_mutation_prob: float=0.2,
                                 new=False):

        # 
        s = np.random.uniform(0,1)

        # the speed is too slow

        if s<apply_mutation_prob:
            if not new:
                mask = self.build_mask(change_percentage)
                chosen = self.dataset[:,mask]
                mutations = np.apply_along_axis(random_substitution, axis=0, arr=chosen)
                self.cell_profile[mask] = mutations[0]
            else:
                mask = self.build_mask(change_percentage)
                cell_random = np.random.randint(self.cell_num, size=int(self.gene_num * change_percentage))
                chosen = self.dataset[cell_random, mask]
                self.cell_profile[mask] = chosen


    def ToTensor(self):
        self.cell_profile = torch.from_numpy(self.cell_profile)


def random_substitution(x):
    random_cell = np.random.randint(x.shape)
    return x[random_cell]
