import argparse
import math
import os
import random
import shutil
import time
import warnings
import numpy as np
import faiss

import torch
import torch.nn as nn
import torch.nn.parallel
import torch.backends.cudnn as cudnn
import torch.optim
import torch.utils.data
import torch.utils.data.distributed
import torchvision.models as models

import pcl.loader
import pcl.builder

from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score
from sklearn.cluster import KMeans

from anndata import AnnData
import scanpy as sc
import pandas as pd

import metrics as M


model_names = sorted(name for name in models.__dict__
    if name.islower() and not name.startswith("__")
    and callable(models.__dict__[name]))

parser = argparse.ArgumentParser(description='PyTorch scRNA-seq CLR Training')

# 1.input h5ad data
parser.add_argument('--input_h5ad_path', type=str, default= "",
                    help='path to input_h5ad')

# 2.hyper-parameters
parser.add_argument('-j', '--workers', default=1, type=int, metavar='N',
                    help='number of data loading workers (default: 32)')

parser.add_argument('--epochs', default=100, type=int, metavar='N',
                    help='number of total epochs to run')

parser.add_argument('--start_epoch', default=0, type=int, metavar='N',
                    help='manual epoch number (useful on restarts)')

parser.add_argument('-b', '--batch_size', default=512, type=int,
                    metavar='N',
                    help='mini-batch size (default: 256), this is the total '
                         'batch size of all GPUs on the current node when '
                         'using Data Parallel or Distributed Data Parallel')

parser.add_argument('--lr', '--learning_rate', default=5e-3, type=float,
                    metavar='LR', help='initial learning rate', dest='lr')

parser.add_argument('--momentum', default=0.9, type=float, metavar='M',
                    help='momentum of SGD solver')

parser.add_argument('--wd', '--weight_decay', default=1e-6, type=float,
                    metavar='W', help='weight decay (default: 1e-4)',
                    dest='weight_decay')

parser.add_argument('--resume', default='', type=str, metavar='PATH',
                    help='path to latest checkpoint (default: none)')

parser.add_argument('--schedule', default=[100, 120], nargs='*', type=int,
                    help='learning rate schedule (when to drop lr by 10x), if use cos, then it will not be activated')

parser.add_argument('--low_dim', default=128, type=int,
                    help='feature dimension (default: 128)')
parser.add_argument('--pcl_r', default=1024, type=int,
                    help='queue size; number of negative pairs; needs to be smaller than num_cluster (default: 16384)')
parser.add_argument('--moco_m', default=0.999, type=float,
                    help='moco momentum of updating key encoder (default: 0.999)')

parser.add_argument('--temperature', default=0.2, type=float,
                    help='softmax temperature')

parser.add_argument('--cos', action='store_true',
                    help='use cosine lr schedule')

parser.add_argument('--warmup_epoch', default=5, type=int,
                    help='number of warm-up epochs to only train with InfoNCE loss')

# augmentation prob
parser.add_argument("--aug_prob", type=float, default=0.5,
                    help="The prob of doing augmentation")

# cluster
parser.add_argument('--cluster_name', default='kmeans', type=str,
                    help='name of clustering method', dest="cluster_name")

parser.add_argument('--num_cluster', default='-1', type=str,
                    help='number of clusters', dest="num_cluster")

# random
parser.add_argument('--seed', default=0, type=int,
                    help='seed for initializing training. ')

# gpu
parser.add_argument('--gpu', default=0, type=int,   #None
                    help='GPU id to use.')

# logs and savings
parser.add_argument('-e', '--eval_freq', default=10, type=int,
                    metavar='N', help='Save frequency (default: 10)',
                    dest='eval_freq')

parser.add_argument('-l', '--log_freq', default=10, type=int,
                    metavar='N', help='print frequency (default: 10)')

parser.add_argument('--exp_dir', default='./experiment_pcl', type=str,
                    help='experiment directory')
# CJY metric
parser.add_argument('--save_dir', default='./result', type=str,
                    help='result saving directory')


def main():
    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)
        torch.manual_seed(args.seed)
        cudnn.deterministic = True
        warnings.warn('You have chosen to seed training. '
                      'This will turn on the CUDNN deterministic setting, '
                      'which can slow down your training considerably! '
                      'You may see unexpected behavior when restarting '
                      'from checkpoints.')

    if args.gpu is not None:
        warnings.warn('You have chosen a specific GPU. This will completely '
                      'disable data parallelism.')
    
    args.num_cluster = args.num_cluster.split(',')

    if not os.path.exists(args.exp_dir):
        os.mkdir(args.exp_dir)

    main_worker(args)


def main_worker(args):
    print(args)

    # 1. Build Dataloader

    # Load h5ad data
    input_h5ad_path = args.input_h5ad_path
    processed_adata = sc.read_h5ad(input_h5ad_path)

    # find dataset name
    pre_path, filename = os.path.split(input_h5ad_path)
    dataset_name, ext = os.path.splitext(filename)
    # for batch effect dataset3
    if dataset_name == "counts":
        dataset_name = pre_path.split("/")[-1]
    if dataset_name == "":
        dataset_name = "unknown"

    # Define Transformation
    args_transformation = {
        # crop
        # without resize, it's better to remove crop
        
        # mask
        'mask_percentage': 0.2,
        'apply_mask_prob': args.aug_prob,
        
        # (Add) gaussian noise
        'noise_percentage': 0.8,
        'sigma': 0.2,
        'apply_noise_prob': args.aug_prob,

        # inner swap
        'swap_percentage': 0.1,
        'apply_swap_prob': args.aug_prob,
        
        # cross over with 1
        'cross_percentage': 0.25,
        'apply_cross_prob': args.aug_prob,
        
        # cross over with many
        'change_percentage': 0.25,
        'apply_mutation_prob': args.aug_prob
    }

    train_dataset = pcl.loader.scRNAMatrixInstance(
        adata=processed_adata,
        transform=True,
        args_transformation=args_transformation
        )
    eval_dataset = pcl.loader.scRNAMatrixInstance(
        adata=processed_adata,
        transform=False
        )

    if train_dataset.num_cells < 512:
        args.batch_size = train_dataset.num_cells
        args.pcl_r = train_dataset.num_cells

    train_sampler = None
    eval_sampler = None
    train_loader = torch.utils.data.DataLoader(
        train_dataset, batch_size=args.batch_size, shuffle=(train_sampler is None),
        num_workers=args.workers, pin_memory=True, sampler=train_sampler, drop_last=True)

    # dataloader for center-cropped images, use larger batch size to increase speed
    eval_loader = torch.utils.data.DataLoader(
        eval_dataset, batch_size=args.batch_size * 5, shuffle=False,
        sampler=eval_sampler, num_workers=args.workers, pin_memory=True)
    
    args.num_cluster = len(train_dataset.unique_label) 

    # 2. Create Model
    print("=> creating model 'MLP'")
    model = pcl.builder.MoCo(
        pcl.builder.MLPEncoder,
        int(train_dataset.num_genes),
        args.low_dim, args.pcl_r, args.moco_m, args.temperature)
    print(model)

    if args.gpu is None:
        raise Exception("Should specify GPU id for training with --gpu".format(args.gpu))
    else:
        print("Use GPU: {} for training".format(args.gpu))

    cudnn.benchmark = True
    torch.cuda.set_device(args.gpu)
    model = model.cuda(args.gpu)
       
    # define loss function (criterion) and optimizer
    criterion = nn.CrossEntropyLoss()   #.cuda(args.gpu)

    optimizer = torch.optim.SGD(model.parameters(), args.lr,
                                momentum=args.momentum,
                                weight_decay=args.weight_decay)

    # optionally resume from a checkpoint
    if args.resume:
        if os.path.isfile(args.resume):
            print("=> loading checkpoint '{}'".format(args.resume))
            if args.gpu is None:
                checkpoint = torch.load(args.resume)
            else:
                # Map model to be loaded to specified single gpu.
                loc = 'cuda:{}'.format(args.gpu)
                checkpoint = torch.load(args.resume, map_location=loc)
            args.start_epoch = checkpoint['epoch']
            model.load_state_dict(checkpoint['state_dict'])
            optimizer.load_state_dict(checkpoint['optimizer'])
            print("=> loaded checkpoint '{}' (epoch {})"
                  .format(args.resume, checkpoint['epoch']))
        else:
            print("=> no checkpoint found at '{}'".format(args.resume))

    # 2. Train Encoder
    # CJY
    # train the model
    for epoch in range(args.start_epoch, args.epochs):

        """
        cluster_result = None
        if epoch >= args.warmup_epoch and epoch % args.save_freq == 0:
            features, labels = inference(eval_loader, model, args)
            if epoch == 10:
                label_decoded = [train_dataset.label_decoder[i] for i in labels]
                df = pd.DataFrame(label_decoded, columns=['x'])
                df.to_csv(os.path.join(args.exp_dir, f'pd_labels_{epoch}.csv'))

            if epoch % args.save_freq==0:# and epoch<300:

                features[np.linalg.norm(features,axis=1)>1.5] /= 2 #account for the few samples that are computed twice  
                
                df = pd.DataFrame(features)
                df.to_csv(os.path.join(args.exp_dir, f'features_{epoch}.csv'))

                #cluster_result = run_kmeans(features,args)  #run kmeans clustering on master node
                # save the clustering result
                #torch.save(cluster_result,os.path.join(args.exp_dir, f'clusters_{epoch}_{args.num_cluster}'))  
                #ari = run_ARI(labels, cluster_result)
                sklearn_kmeans_metrics = run_sklearn_kmeans(features, labels, int(args.num_cluster))

                tmp_ari = sklearn_kmeans_metrics[0]
                tmp_nmi = sklearn_kmeans_metrics[1]
                tmp_pd_label = sklearn_kmeans_metrics[-1]

                if tmp_ari > best_ari:
                    best_ari = tmp_ari
                    best_nmi = tmp_nmi
                    best_features = features
                    best_pd_labels = tmp_pd_label
                    save_labels = labels
            
                # use the leiden algorithm for clustering
                #leiden_ari = run_scanpy_leiden(features, labels)

                with open(os.path.join(args.exp_dir, f'result.txt'), "a") as f:
                    f.writelines(f"{epoch}\t" + '\t'.join((str(elem) for elem in sklearn_kmeans_metrics)) +"\t" + str(acc) + "\n")
                    
        """

        adjust_learning_rate(optimizer, epoch, args)

        # train for one epoch
        train_unsupervised_metrics = train(train_loader, model, criterion, optimizer, epoch, args)

        # training log & unsupervised metrics
        if epoch % args.log_freq == 0 or epoch == args.epochs - 1:
            if epoch == 0:
                with open(os.path.join(args.exp_dir, 'log_{}.txt'.format(dataset_name)), "w") as f:
                    f.writelines(f"epoch\t" + '\t'.join((str(key) for key in train_unsupervised_metrics.keys())) + "\n")
                    f.writelines(f"{epoch}\t" + '\t'.join((str(train_unsupervised_metrics[key]) for key in train_unsupervised_metrics.keys())) + "\n")
            else:
                with open(os.path.join(args.exp_dir, 'log_{}.txt'.format(dataset_name)), "a") as f:
                    f.writelines(f"{epoch}\t" + '\t'.join((str(train_unsupervised_metrics[key]) for key in train_unsupervised_metrics.keys())) + "\n")

        # inference log & supervised metrics
        if epoch % args.eval_freq == 0 or epoch == args.epochs - 1:
            embeddings, gt_labels = inference(eval_loader, model)
            # gt_label exists and metric can be computed
            if gt_labels is not None:
                # perform kmeans
                if args.cluster_name == "kmeans":
                    num_cluster = len(train_dataset.unique_label) if args.num_cluster == -1 else  args.num_cluster
                    pd_labels = KMeans(n_clusters=num_cluster, random_state=args.seed).fit(embeddings).labels_
                    # compute metrics
                    eval_supervised_metrics = M.compute_metrics(gt_labels, pd_labels)
                    print("{}\t {}\n".format(epoch, eval_supervised_metrics))

                    with open(os.path.join(args.exp_dir, 'log_{}.txt'.format(dataset_name)), "a") as f:
                        f.writelines("Eval-epoch-{}\t".format(epoch) + eval_supervised_metrics + "\n")


    # 3. Final Savings
    # save pre-name of dataset
    save_path = args.save_dir

    # save feature & labels
    best_features = embeddings
    best_pd_labels = pd_labels
    best_metrics = eval_supervised_metrics
    gt_labels = gt_labels

    save_path = os.path.join(save_path, "CLEAR")
    if os.path.exists(save_path)!=True:
        os.makedirs(save_path)

    np.savetxt(os.path.join(save_path, "feature_CLEAR_{}.csv".format(dataset_name)), best_features, delimiter=',')
    label_decoded = [train_dataset.label_decoder[i] for i in gt_labels]

    save_labels_df = pd.DataFrame(label_decoded, columns=['x'])
    save_labels_df.to_csv(os.path.join(save_path, "gt_label_CLEAR_{}.csv".format(dataset_name)))

    pd_labels_df = pd.DataFrame(best_pd_labels, columns=['kmeans'])
    pd_labels_df.to_csv(os.path.join(save_path, "pd_label_CLEAR_{}.csv".format(dataset_name)))

    # write metrics into txt
    txt_path = os.path.join(save_path, "metric_CLEAR.txt")
    f = open(txt_path, "a")
    record_string = dataset_name
    for key in best_metrics.keys():
        record_string += " {}".format(best_metrics[key])
    record_string += "\n"
    f.write(record_string)
    f.close()


def train(train_loader, model, criterion, optimizer, epoch, args):
    batch_time = AverageMeter('Time', ':6.3f')
    data_time = AverageMeter('Data', ':6.3f')
    losses = AverageMeter('Loss', ':.4e')
    acc_inst = AverageMeter('Acc@Inst', ':6.2f')   
    
    progress = ProgressMeter(
        len(train_loader),
        [batch_time, data_time, losses, acc_inst],
        prefix="Epoch: [{}]".format(epoch))

    # switch to train mode
    model.train()

    end = time.time()
    for i, (images, index, label) in enumerate(train_loader):
        # measure data loading time
        data_time.update(time.time() - end)

        #import pdb; pdb.set_trace()

        if args.gpu is not None:
            images[0] = images[0].cuda(args.gpu, non_blocking=True)
            images[1] = images[1].cuda(args.gpu, non_blocking=True)
                
        # compute output
        output, target, output_proto, target_proto = model(im_q=images[0], im_k=images[1], cluster_result=None, index=index)
        
        # InfoNCE loss
        loss = criterion(output, target)  

        losses.update(loss.item(), images[0].size(0))
        acc = accuracy(output, target)[0] 
        acc_inst.update(acc[0], images[0].size(0))

        # compute gradient and do SGD step
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        # measure elapsed time 
        batch_time.update(time.time() - end)
        end = time.time()

    progress.display(i+1)

    unsupervised_metrics = {"accuracy": acc_inst.avg.item(), "loss": losses.avg}

    return unsupervised_metrics
            
def inference(eval_loader, model):
    print('Inference...')
    model.eval()
    features = []
    labels = []

    for i, (images, index, label) in enumerate(eval_loader):
        images = images.cuda()
        with torch.no_grad():
            feat = model(images, is_eval=True) 
        feat_pred = feat.data.cpu().numpy()
        label_true = label
        features.append(feat_pred)
        labels.append(label_true)
    
    features = np.concatenate(features, axis=0)
    labels = np.concatenate(labels, axis=0)

    return features, labels

    
def run_sklearn_kmeans(features, labels, num_cluster):
    best_ari = 0
    best_nmi = 0
    best_seed = -1
    for random_seed in range(5):
        kmeans = KMeans(n_clusters=num_cluster, random_state=random_seed).fit(features)
        ari = adjusted_rand_score(labels_true=labels, labels_pred=kmeans.labels_)
        nmi = normalized_mutual_info_score(labels_true=labels, labels_pred=kmeans.labels_)
        if ari > best_ari:
            best_ari = ari
            best_seed = random_seed
        if nmi > best_nmi:
            best_nmi = nmi

    silhou = silhouette_score(features, labels=labels)

    print(f"The sklearn ARI is {best_ari} at sklearn seed {best_seed}\n")
    print(f"The sklearn NMI is {best_nmi} \n")
    print(f"The silhouette coefficient is {silhou}\n")
    
    return [best_ari, best_nmi, silhou, kmeans.labels_]


def run_scanpy_leiden(features, labels):
    adata = AnnData(X=None)
    for i,z in enumerate(features.T):
        adata.obs[f'Z_{i}'] = z
    adata.obsm["X_scVI"] = features
    sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=20)
    #sc.tl.umap(adata, min_dist=0.3)
    sc.tl.leiden(adata, key_added="leiden_scVI", resolution=0.8)


def save_checkpoint(state, is_best, filename='checkpoint.pth.tar'):
    torch.save(state, filename)
    if is_best:
        shutil.copyfile(filename, 'model_best.pth.tar')


class AverageMeter(object):
    """Computes and stores the average and current value"""
    def __init__(self, name, fmt=':f'):
        self.name = name
        self.fmt = fmt
        self.reset()

    def reset(self):
        self.val = 0
        self.avg = 0
        self.sum = 0
        self.count = 0

    def update(self, val, n=1):
        self.val = val
        self.sum += val * n
        self.count += n
        self.avg = self.sum / self.count

    def __str__(self):
        fmtstr = '{name} {val' + self.fmt + '} ({avg' + self.fmt + '})'
        return fmtstr.format(**self.__dict__)


class ProgressMeter(object):
    def __init__(self, num_batches, meters, prefix=""):
        self.batch_fmtstr = self._get_batch_fmtstr(num_batches)
        self.meters = meters
        self.prefix = prefix

    def display(self, batch):
        entries = [self.prefix + self.batch_fmtstr.format(batch)]
        entries += [str(meter) for meter in self.meters]
        print('\t'.join(entries))

    def _get_batch_fmtstr(self, num_batches):
        num_digits = len(str(num_batches // 1))
        fmt = '{:' + str(num_digits) + 'd}'
        return '[' + fmt + '/' + fmt.format(num_batches) + ']'


def adjust_learning_rate(optimizer, epoch, args):
    """Decay the learning rate based on schedule"""
    lr = args.lr
    if args.cos:  # cosine lr schedule
        lr *= 0.5 * (1. + math.cos(math.pi * epoch / args.epochs))
    else:
        # stepwise lr schedule
        for milestone in args.schedule:
            lr *= 0.1 if epoch >= milestone else 1.
    for param_group in optimizer.param_groups:
        param_group['lr'] = lr


def accuracy(output, target, topk=(1,)):
    """Computes the accuracy over the k top predictions for the specified values of k"""
    with torch.no_grad():
        maxk = max(topk)
        batch_size = target.size(0)

        _, pred = output.topk(maxk, 1, True, True)
        pred = pred.t()
        correct = pred.eq(target.view(1, -1).expand_as(pred))

        res = []
        for k in topk:
            correct_k = correct[:k].view(-1).float().sum(0, keepdim=True)
            res.append(correct_k.mul_(100.0 / batch_size))
        return res



if __name__ == '__main__':
    main()
