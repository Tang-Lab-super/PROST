from scipy.sparse import issparse
import random
import numpy as np
import scipy.sparse as sp
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import scanpy as sc

from . utils import setup_seed
from . model import *
from . calculate_PI import *



def run_prost_clust(adata, SEED, n_clusters=None, platform=None, k_neighbors = 7, 
                    min_distance = 50, init="leiden", res = 0.2, key_added = "PROST", 
                    tol=5e-3, gnnlayers = 2, num_pcs = 50, lr = 0.1, dropout=0.05, 
                    leaky_alpha = 0.15, max_epochs = 500, laplacin_filter=True,
                    post_processing = False, pp_run_times = 3):
    
    if init == "louvain" or init == "leiden":
        assert res is not None,"Please set resolution for {} algorithm use ' res = '".format(init)

    else:
        assert n_clusters is not None,"Please set target cluster number for {} algorithm use ' n_clusters = '".format(init)
    
    assert platform is not None,"Please set \" platform='visium' \" or \" platform='slide-seq' \" for run_prost_clust()" 
    setup_seed(SEED)
    
    #--------------------------------------------------------------------------
    print("\nCalculating adjacency matrix ...")
    if platform=="visium":
        adj = get_adj(adata, mode = 'neighbour', n_neighbors = k_neighbors)
    else:
        adj = get_adj(adata, mode = 'distance', min_distance = min_distance)

    #--------------------------------------------------------------------------
    pca = PCA(n_components=num_pcs)
    print("\nRunning PCA ...")
    if sp.issparse(adata.X):
        pca.fit(adata.X.A)
        embed=pca.transform(adata.X.A)
    else:
        pca.fit(adata.X)
        embed=pca.transform(adata.X)
        
    #--------------------------------------------------------------------------
    if laplacin_filter:               
        # embed = sp.csr_matrix(embed).toarray()      
        print('Laplacian Smoothing ...')
        adj_norm = preprocess_graph(adj, gnnlayers, norm='sym', renorm=False) # graph laplacin filter
        for a in adj_norm:
            embed = a.dot(embed)
           
    #--------------------------------------------------------------------------
    Cluster = PROST_cluster()
    Cluster.train(embed, adj,
                    n_pcs = num_pcs, 
                    n_clusters = n_clusters, 
                    init=init, 
                    res=res, 
                    tol=tol, 
                    lr=lr,
                    dropout=dropout,
                    leaky_alpha=leaky_alpha,
                    max_epochs = max_epochs)  
    class_pred, prop, embed = Cluster.predict()
    
    #--------------------------------------------------------------------------
    print("Clustering completed !!")
    adata.obsm[key_added] = embed
    adata.obs["clustering"] = class_pred
    adata.obs["clustering"] = adata.obs["clustering"].astype('category')
    if post_processing:
        cluster_post_process(adata, platform, 
                             k_neighbors = int(k_neighbors*3/2)+1, 
                             min_distance = min_distance*3/2, 
                             key_added = "pp_clustering", 
                             p = 0.5, 
                             run_times = pp_run_times) 
    return adata
            


def cal_prost_index(adata, connect_kernel_size=5, neighbors=8,del_rate=0.01, platform = "visium"):
    adata = minmax_scaler(adata)
    adata = gau_filter(adata, platform)
    adata = get_binary(adata, platform, method = "iterative")
    adata = get_sub(adata, connect_kernel_size, neighbors, platform, del_rate)
    adata = cal_PI(adata,platform)
    print("\nPROST Index calculation completed !!")    
    return adata



