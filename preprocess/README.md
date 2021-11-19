About the scTransform

The python code scTransform.py performs a constraint negative bionomial regression on the scRNA-seq datasets, returns the regressed value for each gene inside the cell. The original paper is here: Hafemeister, C., Satija, R. Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. Genome Biol 20, 296 (2019). https://doi.org/10.1186/s13059-019-1874-1

In our cases, we do augmentations with the regressed values, so during preprocessing, we create 2 AnnData, one for the original data, while the other is for the sctransform version. Both 2 AnnData are required for performing contrastive learning. You can create the sctransformed data by add --scTransform flag in the preprocess step.
