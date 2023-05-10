import pandas as pd
import numpy as np
import os
import math
import random
import numba
import torch
import scipy
import scanpy as sc
from sklearn.neighbors import kneighbors_graph
from sklearn.neighbors import NearestNeighbors
from sklearn import metrics
import scipy.sparse as sp
# from scipy.sparse import issparse
from scipy.optimize import least_squares
from scipy.ndimage import gaussian_filter
from scipy.interpolate import griddata
# from pykrige.ok import OrdinaryKriging
from tqdm import trange
import matplotlib.pyplot as plt



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


def get_adj(adata, mode = "neighbour", n_neighbors = 7, min_distance = 150, self_loop = True):
    spatial = adata.obsm["spatial"]   
    if mode == 'distance':
        assert min_distance is not None,"Please set 'min_diatance' for get_adj()"
        adj = metrics.pairwise_distances(spatial, metric='euclidean')
        adj[adj > min_distance] = 0
        if self_loop:
            adj += np.eye(adj.shape[0])        
        # k = sum(adj>0)  
        return np.int64(adj>0)
    
    elif mode == 'neighbour':
        assert n_neighbors is not None,"Please set 'n_neighbors' for get_adj()"
        adj = kneighbors_graph(spatial, n_neighbors = n_neighbors, include_self = self_loop)
        return adj.toarray()
        

def refine_clusters(result, adj, p):
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



def make_image(genecount, locates, platform = "visium", get_image_idx = False, grid_size = 20,interpolation_method='linear'): # 1d ==> 2d
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

        # 规定插值后尺寸大小
        xloc_new = np.round(locates[:, 0]).astype(int)
        maxx_new, minx_new = np.max(xloc_new), np.min(xloc_new)
        
        yloc_new = np.round(locates[:, 1]).astype(int)
        maxy_new, miny_new = np.max(yloc_new), np.min(yloc_new)

        locates_new = np.array([xloc_new,yloc_new]).T
        #Interpolation
        grid_x, grid_y = np.mgrid[minx_new: maxx_new+1: grid_size, miny_new: maxy_new+1: grid_size]       
        grid_z = griddata(locates, genecount, (grid_x,grid_y), method = interpolation_method) #'nearest''linear''cubic'
        # plt.imshow(grid_z.T)
        # new_shape = grid_z.shape
        return grid_z, grid_z.shape
        

@numba.jit
def gene_img_flatten(I, image_idx_1d): # 2d ==> 1d
    I_1d = I.T.flatten()
    output = np.zeros(image_idx_1d.shape)
    for ii in range(len(image_idx_1d)):
        idx = image_idx_1d[ii]
        output[ii] = I_1d[idx - 1]
    return output


def gau_filter_for_single_gene(gene_data, locates, platform="visium", image_idx_1d=None):
    if platform=="visium":
        I,_ = make_image(gene_data, locates, platform)  
        I = gaussian_filter(I, sigma = 1, truncate = 2)
        output = gene_img_flatten(I, image_idx_1d)
    #--------------------------------------------------------------------------
    else:
        I,_ = make_image(gene_data, locates, platform) 
        # plt.imshow(I.T)
        I = gaussian_filter(I, sigma = 1, truncate = 2)
        # plt.imshow(I.T)
        output = I.flatten()
    return output


def pre_process(adata, percentage = 0.1, min_total_counts = 0, threshold_filter = True, var_stabilization = True):    
    if sp.issparse(adata.X):
        rawcount = adata.X.A.T
    else:
        rawcount = adata.X.T
        
    if percentage > 0:
        count_sum = np.sum(rawcount > 0, 1) #每个基因表达中，计数值大于0的spot数量
        threshold = int(np.size(rawcount, 1) * percentage)
        gene_use = np.where(count_sum >= threshold)[0]
        
    if threshold_filter:
        print("\nFiltering genes ...")
        rawcount = rawcount[gene_use, :]
        
    # if min_total_counts >= 0:
    #     cell_total_counts = np.sum(rawcount, 0)
    #     selects = np.where(cell_total_counts > min_total_counts)[0]
        
    if var_stabilization:
        print("\nVariance-stabilizing transformation to each gene ...")
        rawcount = var_stabilize(rawcount) 
        
    # rawcount = rawcount[:, selects]
  
    return gene_use, rawcount
      
        
def cluster_post_process(adata, platform, k_neighbors, min_distance, 
                         key_added = "pp_clustering", p = 0.5, run_times = 3):
    
    print("\nPost-processing for clustering result ...")
    clutser_result = adata.obs["clustering"]
    # nonlocal PP_adj
    if platform == "visium":
        PP_adj = get_adj(adata, mode = "neighbour", n_neighbors = k_neighbors)
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


def mclust_R(data, num_cluster, modelNames='EEE', random_seed=818):
    """
    Clustering using the mclust algorithm.
    The parameters are the same as those in the R package mclust.
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


def cal_Moran_I(genes_exp, x, y, k=5):
    XYmap = pd.DataFrame({"x": x, "y":y})
    XYnbrs = NearestNeighbors(n_neighbors=k, algorithm='auto', metric = 'euclidean').fit(XYmap)
    XYdistances, XYindices = XYnbrs.kneighbors(XYmap)
    W = np.zeros((genes_exp.shape[0], genes_exp.shape[0]))
    for i in range(0,genes_exp.shape[0]):
        W[i, XYindices[i, :]] = 1
    for i in range(0,genes_exp.shape[0]):
        W[i,i] = 0

    I = pd.Series(index=genes_exp.columns, dtype="float64")
    for k in genes_exp.columns:
        X_minus_mean = np.array(genes_exp[k] - np.mean(genes_exp[k]))
        X_minus_mean = np.reshape(X_minus_mean, (len(X_minus_mean), 1))
        Nom = np.sum(np.multiply(W, np.matmul(X_minus_mean, X_minus_mean.T)))
        Den = np.sum(np.multiply(X_minus_mean, X_minus_mean))
        I[k] = (len(genes_exp[k]) / np.sum(W)) * (Nom / Den)
    return I


def cal_Geary_C(genes_exp, x, y, k=5):
    XYmap = pd.DataFrame({"x": x, "y":y})
    XYnbrs = NearestNeighbors(n_neighbors=k, algorithm='auto', metric = 'euclidean').fit(XYmap)
    XYdistances, XYindices = XYnbrs.kneighbors(XYmap)
    W = np.zeros((genes_exp.shape[0], genes_exp.shape[0]))
    for i in range(0,genes_exp.shape[0]):
        W[i, XYindices[i,:]] = 1
    for i in range(0,genes_exp.shape[0]):
        W[i, i] = 0

    C = pd.Series(index=genes_exp.columns, dtype="float64")
    for k in genes_exp.columns:
        X = np.array(genes_exp[k])
        X_minus_mean = X - np.mean(X)
        X_minus_mean = np.reshape(X_minus_mean, (len(X_minus_mean),1))
        Xij=np.array([X, ] * X.shape[0]).transpose() - np.array([X, ] * X.shape[0])
        Nom = np.sum(np.multiply(W, np.multiply(Xij, Xij)))
        Den = np.sum(np.multiply(X_minus_mean, X_minus_mean))
        C[k] = (len(genes_exp[k]) / (2 * np.sum(W))) * (Nom / Den)
    return C


def cal_moran_I_and_geary_C_for_PI_SVGs(adata, PI_top_n, save_path=None):
    # top_n = 50
    if sp.issparse(adata.X):
        raw_gene_data = adata.X.A.T
    else:
        raw_gene_data = adata.X.T
    pi = pd.DataFrame(adata.var["PI"].reset_index())
    x,y =  adata.obsm["spatial"][:,0], adata.obsm["spatial"][:,1]
    
    # sort data by PI score
    sorted_list = pi.sort_values(by = "PI", ascending = False).reset_index()
    sorted_list = sorted_list.loc[:PI_top_n-1,:]
    sorted_list.columns = ["index","geneID","PI"]
    sorted_index = sorted_list["index"].values    
    selected_raw = pd.DataFrame(raw_gene_data[sorted_index]) 
     
    score_list = sorted_list[["geneID","PI"]]
    score_list["Moran_I"] = score_list["Geary_C"] = np.nan
       
    for i in trange(selected_raw.shape[0]):
        gene_exp =  pd.DataFrame(selected_raw.loc[i,:])

        score_list.loc[i,["Moran_I"]] = cal_Moran_I(gene_exp, x, y).values
        score_list.loc[i,["Geary_C"]] = cal_Geary_C(gene_exp, x, y).values
    
    moran_avg = score_list["Moran_I"].mean()
    geary_avg = score_list["Geary_C"].mean()   
    moran_mid = score_list["Moran_I"].median()
    geary_mid = score_list["Geary_C"].median()
    
    score_list.to_csv("{}/Metrics_of_SVGs.csv".format(save_path))
    print("\nAverage Moran'I of SVGs detected by PI =",moran_avg,
          "\nMedian Moran'I of SVGs detected by PI =",moran_mid,
          "\nAverage Geary'C of SVGs detected by PI =",geary_avg,
          "\nMedian Geary'C of SVGs detected by PI =",geary_mid)
    return score_list
   
    
def get_threshold(adata):
    PI = adata.var["PI"].values
    threshold = np.zeros(3)
    # Probability density distribution of PI value
    kde = scipy.stats.gaussian_kde(PI)
    xs = np.linspace(0, 1, num=1000)
    f = kde(xs)
    x1_idx = np.where(f == max(f))[0].astype(int)
    threshold[0] = xs[x1_idx]
    # First derivative of probability density distribution
    dy = np.diff(f) / np.diff(xs)
    x2_idx1 = np.where(dy == max(dy))[0].astype(int)
    threshold[1] = xs[x2_idx1]
    x3_idx2 = np.where(dy == min(dy))[0].astype(int)
    threshold[2] = xs[x3_idx2]
    return threshold


def preprocess_graph(adj, layer, norm='sym', renorm=True, k = 2/3):
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


def laplacian(adj):
    rowsum = np.array(adj.sum(1))
    degree_mat = sp.diags(rowsum.flatten())
    lap = degree_mat - adj
    return torch.FloatTensor(lap.toarray())


def feature_selection(adata, selected_gene_name=None, by = 'manual', n_top_genes = 3000):
    if by == "prost":
        try:
            pi_score = adata.var["PI"]
        except:
            raise KeyError("Can not find key 'PI' in 'adata.var', please run 'cal_prost_index()' first!")
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
            raise ValueError("Please set 'selected_gene_name' as type [pandas.Series]!")
        selected_genename = pd.Series([i.upper() for i in list(selected_gene_name)])
     
    #
    assert isinstance(selected_genename, pd.Series),"Please input the selected genename as type [pandas.Series]!"
    gene_list = selected_genename.values.tolist()
    gene_list = [i.upper() for i in gene_list]

    adata.var_names = [i.upper() for i in list(adata.var_names)]
    raw_gene_name = adata.var_names.values.tolist()
    for i in raw_gene_name:
        if i in gene_list:
            adata.var.loc[i,['selected']] = 1
        else:
            adata.var.loc[i,['selected']] = 0              
    adata = adata[:, adata.var['selected'] == 1]

    return adata


def cal_metrics_for_DLPFC(labels_pred, labels_true_path, print_result = True):    
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
