import pandas as pd
import numpy as np
import os
import random
import numba
import torch
import scanpy as sc
from sklearn.neighbors import kneighbors_graph
from sklearn.neighbors import NearestNeighbors
from sklearn import metrics
import scipy.sparse as sp
import scipy.stats as stats
from scipy.optimize import least_squares
from scipy.ndimage import gaussian_filter
from scipy.interpolate import griddata
from statsmodels.stats.multitest import fdrcorrection
from tqdm import trange



def setup_seed(seed):
    os.environ['PYTHONHASHSEED'] = str(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.enabled = False


def get_adj(adata, mode = 'neighbour', k_neighbors = 7, min_distance = 150, self_loop = True):
    """
    Calculate adjacency matrix for ST data.
    
    Parameters
    ----------
    mode : str ['neighbour','distance'] (default: 'neighbour')
        The way to define neighbourhood. 
        If `mode='neighbour'`: Calculate adjacency matrix with specified number of nearest neighbors;
        If `mode='distance'`: Calculate adjacency matrix with neighbors within the specified distance.
    k_neighbors : int (default: 7)
        For `mode = 'neighbour'`, set the number of nearest neighbors if `mode='neighbour'`.
    min_distance : int (default: 150)
        For `mode = 'distance'`, set the distance of nearest neighbors if `mode='distance'`.
    self_loop : bool (default: True)
        Whether to add selfloop to the adjacency matrix.
        
    Returns
    -------
    adj : matrix of shape (n_samples, n_samples)
        Adjacency matrix where adj[i, j] is assigned the weight of edge that connects i to j.
    """
    spatial = adata.obsm["spatial"]   
    if mode == 'distance':
        assert min_distance is not None,"Please set `min_diatance` for `get_adj()`"
        adj = metrics.pairwise_distances(spatial, metric='euclidean')
        adj[adj > min_distance] = 0
        if self_loop:
            adj += np.eye(adj.shape[0])  
        adj = np.int64(adj>0)
        return adj
    
    elif mode == 'neighbour':
        assert k_neighbors is not None,"Please set `k_neighbors` for `get_adj()`"
        adj = kneighbors_graph(spatial, n_neighbors = k_neighbors, include_self = self_loop).toarray()
        return adj
        

def var_stabilize(data):
    varx = np.var(data, 1)
    meanx = np.mean(data, 1)
    fun = lambda phi, varx, meanx : meanx + phi * meanx ** 2 - varx
    target_phi = least_squares(fun, x0 = 1, args = (varx, meanx))
    return np.log(data + 1 / (2 * target_phi.x))


def minmax_normalize(data):
    maxdata = np.max(data)
    mindata = np.min(data)
    return (data - mindata)/(maxdata - mindata)


@numba.jit
def get_image_idx_1D(image_idx_2d):
    print("\nCalculating image index 1D:")
    image_idx_1d = np.ones(np.max(image_idx_2d[:])).astype(int)
    for i in trange(1, np.max(image_idx_2d[:])+1):   
        image_idx_1d[i-1] = np.where(image_idx_2d.T.flatten() == i)[0]+1
    return image_idx_1d


def make_image(genecount, locates, platform = "visium", get_image_idx = False, 
               grid_size = 20, interpolation_method='linear'): # 1d ==> 2d
    """
    Convert one-dimensional gene count into two-dimensional interpolated gene image.
    
    Parameters
    ----------
    genecount : pandas.DataFrame
        The matrix of gene count expression. Rows correspond to genes and columns to cells. 
    locates : matrix of shape (n_samples, 2)
        The matrix of gene expression locates. Rows correspond to cells and columns 
        to X-coordinate  and Y-coordinateof the position. 
    platform : str ['visium','Slide-seq','Stereo-seq','osmFISH','SeqFISH'] (default: 'visium')
        Sequencing platforms for generating ST data.
    get_image_idx : bool (default: False)
        If `get_image_idx=True`, calculate `image_idx_1d`. 
        
    grid_size : int (default: 20)
        The size of grid for interpolating irregular spatial gene expression to regular grids.
    interpolation_method : str ['nearest','linear',cubic'] (default: linear)
        The method for interpolating irregular spatial gene expression to regular grids.
        Same as `scipy.interpolate.griddata`
         
    Returns
    -------
    image : ndarray
        2-D gene spatial expression images displayed in a regular pixels.
    image_idx_1d
        If `get_image_idx=True`, which could be input to function `PROST.gene_img_flatten()`.
    """ 
    if platform=="visium":
        xloc = np.round(locates[:, 0]).astype(int)
        maxx = np.max(xloc)
        minx = np.min(xloc)
        yloc = np.round(locates[:, 1]).astype(int)
        maxy = np.max(yloc)
        miny = np.min(yloc)
        
        image = np.zeros((maxy, maxx))    
        image_idx_2d = np.zeros((maxy, maxx)).astype(int)  
        for i in range(len(xloc)):
            temp_y = yloc[i]
            temp_x = xloc[i]
            temp_value = genecount[i]
            image[temp_y - 1, temp_x - 1] = temp_value
            image_idx_2d[temp_y - 1 , temp_x - 1] = i+1
            
        image = np.delete( image, range(miny - 1), 0)
        image = np.delete( image, range(minx - 1), 1)
        image_idx_2d = np.delete(image_idx_2d, range(miny - 1), 0) 
        image_idx_2d = np.delete(image_idx_2d, range(minx - 1), 1)
        image_idx_1d = np.ones(np.max(image_idx_2d[:])).astype(int)
        if get_image_idx:
            image_idx_1d = get_image_idx_1D(image_idx_2d)
                
        return image, image_idx_1d
    #--------------------------------------------------------------------------
    else:
        xloc = locates[:, 0]
        maxx, minx = np.max(xloc), np.min(xloc)

        yloc = locates[:, 1]
        maxy, miny = np.max(yloc), np.min(yloc)

        xloc_new = np.round(locates[:, 0]).astype(int)
        maxx_new, minx_new = np.max(xloc_new), np.min(xloc_new)
        
        yloc_new = np.round(locates[:, 1]).astype(int)
        maxy_new, miny_new = np.max(yloc_new), np.min(yloc_new)

        #Interpolation
        grid_x, grid_y = np.mgrid[minx_new: maxx_new+1: grid_size, miny_new: maxy_new+1: grid_size]       
        image = griddata(locates, genecount, (grid_x,grid_y), method = interpolation_method) #'nearest''linear''cubic'

        return image, image.shape
        

@numba.jit
def gene_img_flatten(I, image_idx_1d): # 2d ==> 1d
    """
    Convert two-dimensional interpolated gene image into one-dimensional gene count.
      
    Parameters
    ----------
    I
        The 2-D gene interpolated image.
    image_idx_1d
        The 2-D index for 1-D gene count. Calculated by function `PROST.make_image()` 
        with setting `get_image_idx = True`

    Returns
    -------
    One-dimensional gene count.
    """ 
    I_1d = I.T.flatten()
    output = np.zeros(image_idx_1d.shape)
    for ii in range(len(image_idx_1d)):
        idx = image_idx_1d[ii]
        output[ii] = I_1d[idx - 1]
    return output


def gau_filter_for_single_gene(gene_data, locates, platform = "visium", image_idx_1d = None):
    """
    Gaussian filter for two-dimensional gene spatial expression images displayed 
    in a regular pixels.
    
    Parameters
    ----------
    gene_data : pandas.DataFrame
        The matrix of gene count expression. Rows correspond to genes and columns to cells. 
    locates : matrix of shape (n_samples, 2)
        The matrix of gene expression locates. Rows correspond to cells and columns 
        to X-coordinate  and Y-coordinateof the position. 
    platform : str ['visium','Slide-seq','Stereo-seq','osmFISH','SeqFISH'] (default: 'visium')
        Sequencing platforms for generating ST data.
    image_idx_1d
        The 2-D index for 1-D gene count. Calculated by function `PROST.make_image()` 
        with setting `get_image_idx = True` 
        
    Returns
    -------
    One-dimensional gene count.
    """ 
    if platform=="visium":
        I,_ = make_image(gene_data, locates, platform)  
        I = gaussian_filter(I, sigma = 1, truncate = 2)
        output = gene_img_flatten(I, image_idx_1d)
    #--------------------------------------------------------------------------
    else:
        I,_ = make_image(gene_data, locates, platform) 
        I = gaussian_filter(I, sigma = 1, truncate = 2)
        output = I.flatten()
    return output


def pre_process(adata, percentage = 0.1, var_stabilization = True):
    """
    Pre-process gene count. 
    
    Parameters
    ----------
    adata : Anndata
        The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond 
        to cells and columns to genes.
    percentage : float (default: 0.1)
        For each gene, count the number of spots (cells) with expression greater 
        than 0, if number is less than a `percentage` of the total spots (cells) 
        number, remove this gene.
    var_stabilization : bool (default: True)
        Var-stabilize transformation.
        
    Returns
    -------
    gene_use
        Index of genes that `percentage` greater than threshold.
    rawcount
        Expression matrix of genes that `percentage` greater than threshold.
    """     
    if sp.issparse(adata.X):
        rawcount = adata.X.A.T
    else:
        rawcount = adata.X.T
        
    if percentage > 0:
        count_sum = np.sum(rawcount > 0, 1) 
        threshold = int(np.size(rawcount, 1) * percentage)
        gene_use = np.where(count_sum >= threshold)[0]
        print("\nFiltering genes ...")
        rawcount = rawcount[gene_use, :]
    else:
        gene_use = np.array(range(len(rawcount)))
          
    if var_stabilization:
        print("\nVariance-stabilizing transformation to each gene ...")
        rawcount = var_stabilize(rawcount) 

    return gene_use, rawcount
      

def refine_clusters(result, adj, p=0.5):
    """
    Reassigning Cluster Labels Using Spatial Domain Information.
    
    Parameters
    ----------
    result
        Clustering result to refine.
    adj
        Adjcency matrix.
    k_neighbors or min_distance
        Different way to calculate adj.
    p : float (default: 0.5)
        Rate of label changes in terms of neighbors
    run_times
        Number of post-process runs. If the label does not change in two consecutive 
        processes, program will also be terminated.
        
    Returns
    -------
    Check post_processed cluster label.
    """
    pred_after = []  
    for i in range(result.shape[0]):
        temp = list(adj[i])  
        temp_list = []
        for index, value in enumerate(temp):
            if value > 0:
                temp_list.append(index) 
        self_pred = result[i]
        neighbour_pred = []      
        for j in temp_list:
            neighbour_pred.append(result[j])
        if (neighbour_pred.count(self_pred) < (len(neighbour_pred))*p) and (neighbour_pred.count(max(set(neighbour_pred), key=neighbour_pred.count))>(len(neighbour_pred))*p):
            pred_after.append(np.argmax(np.bincount(np.array(neighbour_pred))))
        else:
            pred_after.append(self_pred)
    return np.array(pred_after)
      

def cluster_post_process(adata, platform, k_neighbors = None, min_distance = None, 
                         key_added = "pp_clustering", p = 0.5, run_times = 3):
    """
    Post_processing tool for cluster label that integrates neighborhood information.
    
    Parameters
    ----------
    adata : Anndata
        The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond 
        to cells and columns to genes.
    platform : str ['visium','Slide-seq','Stereo-seq','osmFISH','SeqFISH'] (default: 'visium')
        Sequencing platforms for generating ST data.
    k_neighbors : int (default: None)
        Same as `PROST.get_adj()`.
    min_distance : int (default: None)
        Same as `PROST.get_adj()`.
    key_added : str (default: 'pp_clustering')
        `adata.obs` key under which to add the cluster labels.
    p : float (default: 0.5)
        Rate of label changes in terms of neighbors.
    run_times : int (default: 3)
        Number of post-process runs. If the label does not change in two consecutive 
        processes, the run is also terminated.
        
    Returns
    -------
    adata.obs[key_added]
        Array of dim (number of samples) that stores the post-processed cluster 
        label for each cell.
    """
    
    print("\nPost-processing for clustering result ...")
    clutser_result = adata.obs["clustering"]
    # nonlocal PP_adj
    if platform == "visium":
        PP_adj = get_adj(adata, mode = "neighbour", k_neighbors = k_neighbors)
    else:
        PP_adj = get_adj(adata, mode = "distance", min_distance = min_distance)

    result_final = pd.DataFrame(np.zeros(clutser_result.shape[0]))
    i = 1
    while True:        
        clutser_result = refine_clusters(clutser_result, PP_adj, p)
        print("Refining clusters, run times: {}/{}".format(i,run_times))
        result_final.loc[:, i] = clutser_result        
        if result_final.loc[:, i].equals(result_final.loc[:, i-1]) or i == run_times:
            adata.obs[key_added] = np.array(result_final.loc[:, i])
            adata.obs[key_added] = adata.obs[key_added].astype('category')
            return adata
        i += 1


def mclust_R(data, num_cluster, modelNames = 'EEE', random_seed = 818):
    """
    Clustering using the mclust algorithm. 
    Parameters are same as those in the R package mclust.
    """ 
    np.random.seed(random_seed)
    import rpy2.robjects as robjects
    robjects.r.library("mclust")
    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()  
    r_random_seed = robjects.r['set.seed']
    r_random_seed(random_seed)
    rmclust = robjects.r['Mclust']
    res = rmclust(data, num_cluster, modelNames)
    return np.array(res[-2])


def calc_I(w, y):
    """
    Calculate Moran's I.
    
    Parameters
    ----------
    w
        Spatial weights.
    y
        Attribute vector.

    Returns
    -------
    The value of Moran's I
    """ 
    n = len(y)
    z = y - y.mean()
    z = z.reshape(len(z),1)
    zl = w * z
    num = (zl * z.T).sum()
    z2ss = (z * z).sum()
    return n / w.sum() * num / z2ss


def calc_C(w, y):
    """
    Calculate Geary's C.
    
    Parameters
    ----------
    w
        Spatial weights.
    y
        Attribute vector.
    
    Returns
    -------
    The value of Geary's C
    """ 
    n = len(y)
    s0 = w.sum()
    z = y - y.mean()
    z = z.reshape(len(z),1)
    z2ss = (z * z).sum()
    den = z2ss * s0 * 2.0
    num = (w * ( (np.array([y, ] * y.shape[0]).transpose() - np.array([y, ] * y.shape[0]))** 2)).sum()
    return (n - 1) * num / den


def spatial_autocorrelation(adata, k = 10, permutations = None):
    """
    Statistical test of spatial autocorrelation for each gene.
    
    Parameters
    ----------
    adata : Anndata
        The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.
    k : int (default: 10)
        Number of neighbors to define neighborhood.
    permutations : int (default: None)
        Number of random permutations for calculating pseudo p-values. 
        Default is 'none' to skip this step.
        
    Returns
    -------
    adata : Anndata
        adata.var["Moran_I"] : Moran's I
        adata.var["Geary_C"] : Moran's C
        adata.var["p_norm"] : p-value under normality assumption
        adata.var["p_rand"] : p-value under randomization assumption
        adata.var["fdr_norm"] : FDR under normality assumption
        adata.var["fdr_rand"] : FDR under randomization assumption
        
        if set `permutations`:
        adata.var["p_sim"] : p-value based on permutation test
        adata.var["fdr_sim"] : FDR based on permutation test
    """
    if sp.issparse(adata.X):
        genes_exp = adata.X.A
    else:
        genes_exp = adata.X
    spatial = adata.obsm['spatial'] 
    w = kneighbors_graph(spatial, n_neighbors = k, include_self = False).toarray()
    
    s0 = w.sum()
    s02 = s0 * s0
    t = w + w.transpose()
    s1 = np.multiply(t, t).sum()/2.0
    s2 = (np.array(w.sum(1) + w.sum(0).transpose()) ** 2).sum()
    n = len(genes_exp)
    n2 = n * n
    E = -1.0 / (n - 1)
    
    v_num = n2 * s1 - n * s2 + 3 * s02
    v_den = (n - 1) * (n + 1) * s02
    V_norm = v_num / v_den - (1.0 / (n - 1)) ** 2
    
    moranI = np.zeros(genes_exp.shape[1])
    gearyC = np.zeros(genes_exp.shape[1])
    p_norm = np.zeros(genes_exp.shape[1])
    p_rand = np.zeros(genes_exp.shape[1])
    p_sim = np.zeros(genes_exp.shape[1])
    for i in trange(genes_exp.shape[1]):     
        exp = genes_exp[:, i]
        moranI[i] = calc_I(w, exp)
        gearyC[i] = calc_C(w, exp)
        
        z = exp - exp.mean()
        z2 = z ** 2
        z4 = z ** 4
        D = (z4.sum() / n) / ((z2.sum() / n) ** 2)
        A = n * ((n2 - 3 * n + 3) * s1 - n * s2 + 3 * s02)
        B = D * ((n2 - n) * s1 - 2 * n * s2 + 6 * s02)
        C = ((n - 1) * (n - 2) * (n - 3) * s02)
        E_2 = (A - B) / C
        V_rand = E_2 - E * E
        
        z_norm = (moranI[i]-E) / V_norm**(1 / 2.0)
        z_rand = (moranI[i]-E) / V_rand**(1 / 2.0)
        
        p_norm[i] = stats.norm.sf(abs(z_norm))
        p_rand[i] = stats.norm.sf(abs(z_rand))
            
        if permutations:
            sim = [ calc_I(w, np.random.permutation(exp)) for i in range(permutations)]
            sim = np.array(sim)
            larger = np.sum(sim>=moranI[i])

            if (permutations - larger) < larger:
                larger = permutations - larger
            p_sim[i] = larger+1 / permutations+1
            
    _, fdr_norm = fdrcorrection(p_norm, alpha=0.05)
    _, fdr_rand = fdrcorrection(p_rand, alpha=0.05)
    
    adata.var["Moran_I"] = moranI
    adata.var["Geary_C"] = gearyC
    adata.var["p_norm"] = p_norm # p-value under normality assumption
    adata.var["p_rand"] = p_rand # p-value under randomization assumption
    adata.var["fdr_norm"] = fdr_norm
    adata.var["fdr_rand"] = fdr_rand

    if permutations:
        _, fdr_sim = fdrcorrection(p_sim, alpha=0.05)
        adata.var["p_sim"] = p_sim
        adata.var["fdr_sim"] = fdr_sim
    
    return adata


def preprocess_graph(adj, layer = 2, norm = 'sym', renorm = True, k = 2/3):
    """
    Preprocess adj matrix.
    
    Parameters
    ----------
    adata : Anndata
        The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond 
        to cells and columns to genes.
    selected_gene_name : pd.Series (default: None)
        Manually set the genes' name to seclect. Input as type [pandas.Series]
    by : str ["prost", "scanpy"] (default: None)
        Method for feature selection. 
        If `by=="prost"`, feature will be selected by PI;
        If `by=="scanpy"`, feature will be selected by Seurat.
    n_top_genes : int (default: 3000)
        Number of features (spatially variable genes) to select.
        
    Returns
    -------
    adata that include only selected genes.
    """
    adj = sp.coo_matrix(adj)
    ident = sp.eye(adj.shape[0])
    if renorm:
        adj_ = adj + ident
    else:
        adj_ = adj  
    rowsum = np.array(adj_.sum(1)) 
    
    if norm == 'sym':
        degree_mat_inv_sqrt = sp.diags(np.power(rowsum, -0.5).flatten())
        adj_normalized = adj_.dot(degree_mat_inv_sqrt).transpose().dot(degree_mat_inv_sqrt).tocoo()
        laplacian = ident - adj_normalized
        
    elif norm == 'left':
        degree_mat_inv_sqrt = sp.diags(np.power(rowsum, -1.).flatten())
        adj_normalized = degree_mat_inv_sqrt.dot(adj_).tocoo()
        laplacian = ident - adj_normalized  
        
    reg = [k] * layer
    adjs = []
    for i in range(len(reg)):
        adjs.append(ident-(reg[i] * laplacian))
    return adjs


def feature_selection(adata, selected_gene_name = None, by = 'prost', n_top_genes = 3000):
    """
    A feature selection tool for ST data.
    
    Parameters
    ----------
    adata : Anndata
        The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond 
        to cells and columns to genes.
    selected_gene_name : pd.Series (default: None)
        Manually set the genes' name to seclect. Input as type [pandas.Series]
    by : str ["prost", "scanpy"] (default: None)
        Method for feature selection. 
        If `by=="prost"`, feature will be selected by PI;
        If `by=="scanpy"`, feature will be selected by Seurat.
    n_top_genes : int (default: 3000)
        Number of features (spatially variable genes) to select.
        
    Returns
    -------
    adata that include only selected genes.
    """
    if by == "prost":
        try:
            pi_score = adata.var["PI"]
        except:
            raise KeyError("Can not find key 'PI' in 'adata.var', please run `PROST.cal_prost_index()` first!")
        pi_score = adata.var["PI"]
        sorted_score = pi_score.sort_values(ascending = False)
        gene_num = np.sum(sorted_score>0)
        selected_num = np.minimum(gene_num,n_top_genes)
        selected_genename = pd.Series(sorted_score[:selected_num].index)          
    elif by == "scanpy":
        sc.pp.highly_variable_genes(adata, n_top_genes = n_top_genes)
        adata = adata[:, adata.var.highly_variable]
        return adata
    else:
        if selected_gene_name is None:
            raise ValueError("Please set 'selected_gene_name' as type `pandas.Series` !")
        selected_genename = pd.Series([i.upper() for i in list(selected_gene_name)])
     
    #
    assert isinstance(selected_genename, pd.Series),"Please input the `selected_gene_name` as type `pandas.Series`!"
    gene_list = selected_genename.values.tolist()
    gene_list = [i.upper() for i in gene_list]

    adata.var_names = [i.upper() for i in list(adata.var_names)]
    raw_gene_name = adata.var_names.values.tolist()
    for i in raw_gene_name:
        if i in gene_list:
            adata.var.loc[i,['selected']] = True
        else:
            adata.var.loc[i,['selected']] = False  
    adata.var.selected = adata.var.selected.astype('bool')
    adata = adata[:, adata.var.selected]

    return adata


def cal_metrics_for_DLPFC(labels_pred, labels_true_path=None, print_result = True):

    # labels_true processing
    labels_true = pd.read_csv(labels_true_path)
    labels_true['ground_truth'] = labels_true['ground_truth'].str[-1]
    labels_true = labels_true.fillna(8)   
    for i in range(labels_true.shape[0]):
        temp = labels_true['ground_truth'].iloc[i]
        if temp == 'M':
            labels_true['ground_truth'].iloc[i] = 7       
    labels_true = pd.DataFrame(labels_true['ground_truth'], dtype=np.int64).values
    labels_true = labels_true[:,0]    
    #
    ARI = metrics.adjusted_rand_score(labels_true, labels_pred)
    AMI = metrics.adjusted_mutual_info_score(labels_true, labels_pred)
    NMI = metrics.normalized_mutual_info_score(labels_true, labels_pred)
    v_measure_score = metrics.v_measure_score(labels_true, labels_pred)
    silhouette_score = metrics.silhouette_score(np.array(labels_true).reshape(-1, 1), np.array(labels_pred).reshape(-1, 1).ravel())
    if print_result:
        print('\nARI =', ARI, '\nAMI =', AMI, '\nNMI =', NMI, 
              '\nv_measure_score =', v_measure_score, '\nsilhouette_score =',silhouette_score,
              '\n==================================================================')
    return ARI, NMI, silhouette_score


def simulateH5Data(adata, rr = 0.05):
    """
    Get simulated data by droping part of the gene randomly.
    
    Parameters
    ----------
    adata : Anndata
        H5 object. The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.
    rr : float (default: 0,05)
        Dropout rate.
        
    Returns
    -------
    adata with specified dropout rates.
    """
    if rr > 1 or rr < 0:
        print("Warning! Dropout rate is illegal!")
        return 0
    print(f"\ndropout rate = {rr}")
    
    from random import sample
    import copy

    # get expression matrix
    issparse = 0
    if sp.issparse(adata.X):
        data_ori_dense = adata.X.A
        issparse = 1
    else:
        data_ori_dense = adata.X
    
    # sample from non-zero
    flagXY = np.where(data_ori_dense != 0)      
    ncount = len(flagXY[0])

    # sample rr% -> 0.0
    flag = sample(range(ncount), k=int(rr*ncount))
    dropX, dropY = flagXY[0][flag], flagXY[1][flag]

    # update anndata
    data_new = data_ori_dense.copy()
    for dx, dy in zip(dropX, dropY):
        data_new[dx, dy] = 0.0
    reCount = (data_new != 0).sum()
    if issparse:
        data_new = sp.csr_matrix(data_new)

    # new adata, return
    newAdata = copy.deepcopy(adata)
    newAdata.X = data_new

    # Note:not Updata metadata!
    print(f"Done! Remain {reCount}/{ncount}") 
    return newAdata