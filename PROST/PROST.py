import torch
import numpy as np
import scipy.sparse as sp
from sklearn.decomposition import PCA

from . utils import setup_seed, get_adj, preprocess_graph, cluster_post_process
from . model import PROST_NN, PROST_NN_sparse
from . calculate_PI import minmax_scaler, gau_filter, get_binary, get_sub, cal_prost_index



def cal_PI(adata, kernel_size=5, del_rate=0.01, platform = "visium"):
    '''
    Use PI to identify spatially variable genes for ST data.
    
    Parameters
    ----------
    adata : Anndata
        The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    grid_size : int (default: 5)
        Grid size for interpolation. Set this parameter only when `platform != visium`
    kernel_size : int (default: 5)
        Define the size of kernel for implementing morphological closing operation 
        on the binary image. The larger size of kernel could eliminate a larger 
        minutiae section inside the foreground.      
    del_rate : float (default: 0.01)
        For a gene, if the largest foreground is less than `del_rate` of the 
        entire spatial expression, it will be recognized as not having a significant 
        spatial pattern.       
    platform : str ['visium','Slide-seq','Stereo-seq','osmFISH','SeqFISH' or other platform that generate irregular spots] (default: 'visium')
        Sequencing platforms for generating ST data.
        
    Returns
    -------
    adata : Anndata
        adata.obs.PI : PI score for each genes.
        adata.obs.SEP : Seperability score for each genes.
        adata.obs.SIG : Significance score for each genes.
    '''
    
    adata = minmax_scaler(adata)
    adata = gau_filter(adata, platform)
    adata = get_binary(adata, platform, method = "iterative")
    adata = get_sub(adata, kernel_size, platform, del_rate)
    adata = cal_prost_index(adata, platform)
    print("\nPROST Index calculation completed !!")    
    return adata


def run_PNN(adata, SEED, platform=None, init="leiden", n_clusters=5, res=0.2,
                    k_neighbors=7, min_distance=50, key_added="PROST", lap_filter=2, 
                    lr=0.1, tol=5e-3, max_epochs=500, post_processing=False, 
                    pp_run_times=3, cuda=False):
    '''
    Use PNN to identify spatial domains for ST data.
    
    Parameters
    ----------
    adata : Anndata
        The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.
    SEED : int
        Random seed.
    platform : str ['visium','Slide-seq','Stereo-seq','osmFISH','SeqFISH' or other platform that generate irregular spots] (default: 'visium')
        Sequencing platforms for generating ST data.
    init : str ["kmeans","mclust","louvain","leiden"] (default: leiden)
        Methods for initializing cluster centroids.
    n_clusters : int (default: 5)
        If the number of spatial domains is know, set cluster numbers for `init='kmeans'` or `init='mclust'`.
    res : float (default: 0.5)
        If the number of spatial domains is unknown, set resolutions parameter for `init='kmeans'` or `init='mclust'`.
    k_neighbors : int (default: 7)
        For `mode = 'neighbour'`, set the number of nearest neighbors if `mode='neighbour'`.
    min_distance : int (default: 50)
        For `mode = 'distance'`, set the distance of nearest neighbors if `mode='distance'`.   
    key_added : str (default: 'PROST')
        `adata.obsm` key under which to add the embedding representation generated by PROST.
    lap_filter : int (default: 2)
        Number of stacked laplacian filter.
    lr : float (default: 0.1)
        Learning rate.
    tol : float (default: 5e-3)
        Stop criterion. The procedure stops when the change of clustering assignment between two consecutive iterations less than `tol`.
    max_epochs : int (default: 500)
        Number of epoch to train model.
    post_processing : bool (default: False)
        Whether to post-process oringal cluster result.
    pp_run_times : int (default: 3)
        For `post_processing=True`, set the number of post-processing run times.
    cuda : bool (default: False)
        Whether to use cuda acceleration.
    
    Returns
    -------
    adata : Anndata
        adata.obs['clustering'] : Original cluster label for each spot(cell).
        adata.obs['pp_clustering'] : Post-processed cluster label for each spot(cell).
        adata.obsm['PROST'] : Embedding representation generated by PROST.  
    '''

    setup_seed(SEED)
 
    #--------------------------------------------------------------------------
    print("\nCalculating adjacency matrix ...")
    if platform=="visium" or platform=="ST":
        adj = get_adj(adata, mode = 'neighbour', k_neighbors = k_neighbors)
        adj = adj.toarray()
    else:
        adj = get_adj(adata, mode = 'distance', min_distance = min_distance)

    #--------------------------------------------------------------------------
    num_pcs = int(min(50, adata.shape[1]))
    pca = PCA(n_components=num_pcs)
    print("\nRunning PCA ...")
    if sp.issparse(adata.X):
        pca.fit(adata.X.A)
        embed=pca.transform(adata.X.A)
    else:
        pca.fit(adata.X)
        embed=pca.transform(adata.X)
        
    #--------------------------------------------------------------------------
    if lap_filter>0 and lap_filter!=False:
        print('Laplacian Smoothing ...') # Graph laplacin smoothing 
        adj_norm = preprocess_graph(adj, lap_filter, norm='sym', renorm=False)
        for a in adj_norm:
            embed = a.dot(embed)

    #--------------------------------------------------------------------------
    PNN = PROST_NN(embed.shape[1], embed.shape[1], cuda)
    
    PNN.train_(embed, adj, init=init, n_clusters=n_clusters, res=res, tol=tol, 
               lr=lr, max_epochs=max_epochs, seed=SEED)  
    
    embed, prop = PNN.predict(embed, adj)
    print("Clustering completed !!")
    
    #--------------------------------------------------------------------------
    embed = embed.detach().cpu().numpy()
    y_pred = torch.argmax(prop, dim=1).data.cpu().numpy()   
    
    adata.obsm[key_added] = embed
    adata.obs["clustering"] = y_pred
    adata.obs["clustering"] = adata.obs["clustering"].astype('category')

    #--------------------------------------------------------------------------
    if post_processing:
        cluster_post_process(adata, platform, 
                             k_neighbors = int(k_neighbors*3/2)+1, 
                             min_distance = min_distance*3/2, 
                             key_added = "pp_clustering", 
                             run_times = pp_run_times) 
    return adata
            

def run_PNN_sparse(adata, SEED, init="leiden", n_clusters=5, res=0.5,
                    k_neighbors=7, key_added="PROST", lap_filter=2, 
                    lr=0.1, tol=5e-3, max_epochs=500, cuda=False):
    '''
    Sparse version of PNN.
    
    Parameters
    ----------
    adata : Anndata
        The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.
    SEED : int
        Random seed.
    init : str ["kmeans","mclust","louvain","leiden"] (default: leiden)
        Methods for initializing cluster centroids.
    n_clusters : int (default: 5)
        If the number of spatial domains is know, set cluster numbers for `init='kmeans'` or `init='mclust'`.
    res : float (default: 0.5)
        If the number of spatial domains is unknown, set resolutions parameter for `init='kmeans'` or `init='mclust'`.
    k_neighbors : int (default: 7)
        Set the number of nearest neighbors to calculate adjcency matrix.
    key_added : str (default: 'PROST')
        `adata.obsm` key under which to add the embedding representation generated by PROST.
    lap_filter : int (default: 2)
        Number of stacked laplacian filter.
    lr : float (default: 0.1)
        Learning rate.
    tol : float (default: 5e-3)
        Stop criterion. The procedure stops when the change of clustering assignment between two consecutive iterations less than `tol`.
    max_epochs : int (default: 500)
        Number of epoch to train model.
    cuda : bool (default: False)
        Whether to use cuda acceleration.
    
    Returns
    -------
    adata : Anndata
        adata.obs['clustering'] : Original cluster label for each spot(cell).
        adata.obs['pp_clustering'] : Post-processed cluster label for each spot(cell).
        adata.obsm['PROST'] : Embedding representation generated by PROST.  
    '''

    setup_seed(SEED)
 
    #--------------------------------------------------------------------------
    print("\nCalculating adjacency matrix ...")
    adj = get_adj(adata, mode = 'neighbour', k_neighbors = k_neighbors)

    #--------------------------------------------------------------------------
    num_pcs = int(min(50, adata.shape[1]))
    pca = PCA(n_components=num_pcs)
    print("\nRunning PCA ...")
    if sp.issparse(adata.X):
        pca.fit(adata.X.A)
        embed=pca.transform(adata.X.A)
    else:
        pca.fit(adata.X)
        embed=pca.transform(adata.X)
        
    #--------------------------------------------------------------------------
    print('Laplacian Smoothing ...') # Graph laplacin smoothing 
    adj_norm = preprocess_graph(adj, lap_filter, norm='sym', renorm=False) 
    for a in adj_norm:
        embed = a.dot(embed)
           
    #--------------------------------------------------------------------------
    PNN = PROST_NN_sparse(embed.shape[1], embed.shape[1], cuda)
    
    PNN.train_(embed, adj, init=init, n_clusters=n_clusters, res=res, tol=tol, 
               lr=lr, max_epochs=max_epochs, seed=SEED)  
    
    embed, prop = PNN.predict(embed, adj)
    print("Clustering completed !!")
    
    #--------------------------------------------------------------------------
    embed = embed.detach().cpu().numpy()
    y_pred = torch.argmax(prop, dim=1).data.cpu().numpy()   
    
    adata.obsm[key_added] = embed
    adata.obs["clustering"] = y_pred
    adata.obs["clustering"] = adata.obs["clustering"].astype('category')

    return adata
