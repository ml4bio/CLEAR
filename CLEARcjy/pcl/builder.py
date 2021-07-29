import torch
import torch.nn as nn
from random import sample
from torch.nn.modules.linear import Linear
import torchvision.models as models
import math
import torch.nn.functional as F
from torch.nn.parameter import Parameter








def full_block(in_features, out_features, p_drop=0.0):
    return nn.Sequential(
        nn.Linear(in_features, out_features, bias=True),
        #nn.LayerNorm(out_features),
        nn.ReLU(),
        nn.Dropout(p=p_drop),
    )






class MLPEncoder(nn.Module):

    def __init__(self, num_genes=10000, num_hiddens=128, p_drop=0.0):
        super().__init__()
        self.encoder = nn.Sequential(
            full_block(num_genes, 1024, p_drop),
            full_block(1024, num_hiddens, p_drop),
            # add one block for features
        )

    
    def forward(self, x):

        x = self.encoder(x)

        return x


class Classification_model(nn.Module):
    def __init__(self, 
                 num_genes,
                 num_hiddens,
                 p_drop,
                 num_clusters):
        super().__init__()
        self.encoder = MLPEncoder(num_genes, num_hiddens, p_drop)
        self.classification = full_block(num_hiddens, num_clusters)

    def forward(self, x):
        feat = self.encoder(x)
        x = self.classification(feat)
        return feat, x


class MoCo(nn.Module):
    """
    Build a MoCo model with: a query encoder, a key encoder, and a queue
    https://arxiv.org/abs/1911.05722
    """
    def __init__(self, base_encoder, num_genes=10000,  dim=16, r=512, m=0.999, T=0.2):
        """
        dim: feature dimension (default: 16)
        r: queue size; number of negative samples/prototypes (default: 512)
        m: momentum for updating key encoder (default: 0.999)
        T: softmax temperature 
        mlp: whether to use mlp projection
        """
        super(MoCo, self).__init__()

        self.r = r
        self.m = m
        self.T = T

        # create the encoders
        # num_classes is the output fc dimension
        self.encoder_q = base_encoder(num_genes=num_genes, num_hiddens=dim)
        self.encoder_k = base_encoder(num_genes=num_genes, num_hiddens=dim)

        for param_q, param_k in zip(self.encoder_q.parameters(), self.encoder_k.parameters()):
            param_k.data.copy_(param_q.data)  # initialize
            param_k.requires_grad = False  # not update by gradient

        # create the queue
        self.register_buffer("queue", torch.randn(dim, r))
        self.queue = nn.functional.normalize(self.queue, dim=0)

        self.register_buffer("queue_ptr", torch.zeros(1, dtype=torch.long))

    @torch.no_grad()
    def _momentum_update_key_encoder(self):
        """
        Momentum update of the key encoder
        """
        for param_q, param_k in zip(self.encoder_q.parameters(), self.encoder_k.parameters()):
            param_k.data = param_k.data * self.m + param_q.data * (1. - self.m)

    @torch.no_grad()
    def _dequeue_and_enqueue(self, keys):
        # gather keys before updating queue
        # keys = concat_all_gather(keys)

        batch_size = keys.shape[0]

        ptr = int(self.queue_ptr)
        assert self.r % batch_size == 0  # for simplicity

        # replace the keys at ptr (dequeue and enqueue)
        self.queue[:, ptr:ptr + batch_size] = keys.T
        ptr = (ptr + batch_size) % self.r  # move pointer

        self.queue_ptr[0] = ptr

    # @torch.no_grad()
    # def _batch_shuffle_ddp(self, x):
    #     """
    #     Batch shuffle, for making use of BatchNorm.
    #     *** Only support DistributedDataParallel (DDP) model. ***
    #     """
    #     # gather from all gpus
    #     batch_size_this = x.shape[0]
    #     x_gather = concat_all_gather(x)
    #     batch_size_all = x_gather.shape[0]

    #     num_gpus = batch_size_all // batch_size_this

    #     # random shuffle index
    #     idx_shuffle = torch.randperm(batch_size_all).cuda()

    #     # broadcast to all gpus
    #     torch.distributed.broadcast(idx_shuffle, src=0)

    #     # index for restoring
    #     idx_unshuffle = torch.argsort(idx_shuffle)

    #     # shuffled index for this gpu
    #     gpu_idx = torch.distributed.get_rank()
    #     idx_this = idx_shuffle.view(num_gpus, -1)[gpu_idx]

    #     return x_gather[idx_this], idx_unshuffle

    # @torch.no_grad()
    # def _batch_unshuffle_ddp(self, x, idx_unshuffle):
    #     """
    #     Undo batch shuffle.
    #     *** Only support DistributedDataParallel (DDP) model. ***
    #     """
    #     # gather from all gpus
    #     batch_size_this = x.shape[0]
    #     x_gather = concat_all_gather(x)
    #     batch_size_all = x_gather.shape[0]

    #     num_gpus = batch_size_all // batch_size_this

    #     # restored index for this gpu
    #     gpu_idx = torch.distributed.get_rank()
    #     idx_this = idx_unshuffle.view(num_gpus, -1)[gpu_idx]

    #     return x_gather[idx_this]

    def forward(self, im_q, im_k=None, is_eval=False, cluster_result=None, index=None):
        """
        Input:
            im_q: a batch of query images
            im_k: a batch of key images
            is_eval: return momentum embeddings (used for clustering)
            cluster_result: cluster assignments, centroids, and density
            index: indices for training samples
        Output:
            logits, targets, proto_logits, proto_targets
        """
        
        if is_eval:
            k = self.encoder_k(im_q)  
            k = nn.functional.normalize(k, dim=1)            
            return k
        
        # compute key features
        with torch.no_grad():  # no gradient to keys
            self._momentum_update_key_encoder()  # update the key encoder

            # shuffle for making use of BN
            # im_k, idx_unshuffle = self._batch_shuffle_ddp(im_k)

            k = self.encoder_k(im_k)  # keys: NxC
            # print(k.shape)
            k = nn.functional.normalize(k, dim=1)

            # undo shuffle
            #k = self._batch_unshuffle_ddp(k, idx_unshuffle)

        # compute query features
        q = self.encoder_q(im_q)  # queries: NxC
        q = nn.functional.normalize(q, dim=1)
        
        # compute logits
        # Einstein sum is more intuitive
        # positive logits: Nx1
        l_pos = torch.einsum('nc,nc->n', [q, k]).unsqueeze(-1)
        # negative logits: Nxr
        l_neg = torch.einsum('nc,ck->nk', [q, self.queue.clone().detach()])

        # logits: Nx(1+r)
        logits = torch.cat([l_pos, l_neg], dim=1)

        # apply temperature
        logits /= self.T

        # labels: positive key indicators
        labels = torch.zeros(logits.shape[0], dtype=torch.long).cuda()

        # dequeue and enqueue
        self._dequeue_and_enqueue(k)
        

        return logits, labels, None, None



