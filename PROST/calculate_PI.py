import numpy as np
import cv2
import scipy
import scipy.sparse as sp
import scanpy as sc
from tqdm import trange
from skimage.measure import label
from . utils import *



def prepare_for_PI(adata, percentage = 0.1, platform="visium"):    
    selected_gene_idxs, postcount = pre_process(adata, percentage, threshold_filter = True, var_stabilization = False)  
    # _, varcount = pre_process(adata, threshold_filter = False, var_stabilization = True) 
    
    if platform=="visium":
        # _, varcount = pre_process(adata, threshold_filter = False, var_stabilization = True) 
        locates = adata.obs[["array_row","array_col"]].values
        if np.min(locates) == 0:
            locates += 1
        _, image_idx = make_image(postcount[0], locates, platform, get_image_idx = True)

        # selected_var = varcount[selected_gene_idxs, :]
        adata = adata[:, selected_gene_idxs]
        sc.pp.filter_genes(adata, min_cells=3)
        adata.obs['image_idx_1d'] = image_idx
        # adata.uns['selected_var'] = selected_var
    else:
        locates = adata.obsm["spatial"].astype(float)
        _, shape = make_image(postcount[0], locates, platform)
        adata = adata[:, selected_gene_idxs]
        sc.pp.filter_genes(adata, min_cells=50)
        adata.uns['shape'] = shape

    adata.uns['locates'] = locates

    return adata


def minmax_scaler(adata):
    if sp.issparse(adata.X):
        data = adata.X.A.T
    else:
        data = adata.X.T        
    print('\nNormalization to each gene:')
    nor_counts = np.zeros(data.shape)
    for i in trange(len(data)):
        nor_counts[i, :] = minmax_normalize(data[i])    
    adata.uns['nor_counts'] = nor_counts  
    return adata


def gau_filter(adata,platform="visium"):
    gene_data = adata.uns['nor_counts']
    locates = adata.uns['locates']
    print('\nGaussian filtering for each gene:')
    
    #--------------------------------------------------------------------------
    if platform=="visium":
        image_idx_1d = adata.obs['image_idx_1d'].astype(int)
        gau_fea = np.zeros(gene_data.shape)
        for i in trange(len(gene_data)):
            gau_fea[i, :] = gau_filter_for_single_gene(gene_data[i], locates, platform, image_idx_1d)
            
    #--------------------------------------------------------------------------
    else:
        gau_fea = np.zeros((gene_data.shape[0],adata.uns['shape'][0]*adata.uns['shape'][1]))
        for i in trange(len(gene_data)):
            gau_fea[i, :] = gau_filter_for_single_gene(gene_data[i], locates, platform)
    # print('Saveing...')       
    adata.uns['gau_fea'] = gau_fea     
    return adata


def get_binary(adata, platform="visium", method = "iterative"):
    gene_data = adata.uns['gau_fea']
    if sp.issparse(adata.X):
        raw_gene_data = adata.X.A.T
    else:
        raw_gene_data = adata.X.T
    locates = adata.uns['locates']
    
    print('\nBinary segmentation for each gene:')
    
    #--------------------------------------------------------------------------
    if platform=="visium":
        image_idx_1d = adata.obs['image_idx_1d']
        output = np.zeros(gene_data.shape)  
        for i in trange(len(gene_data)):
            fig1 = gene_data[i, :]
            # r1 = raw_gene_data[i, :]
            Im, _ = make_image(fig1, locates, platform)
            if method == "iterative":    
                    m, n = Im.shape       
                    zd = float(np.max(Im))
                    zx = float(np.min(Im))
                    Th = float((zd+zx))/2
                    while True:
                        S0 = 0.0; n0 = 0.0; S1 = 0.0; n1 = 0.0
                        for p in range(m):
                            for q in range(n):
                                if float(Im[p,q]) >= Th:
                                    S1 += float(Im[p,q])
                                    n1 += 1
                                else:
                                    S0 += float(Im[p,q])
                                    n0 += 1
                        T0 = S0/n0; T1 = S1/n1
                        if abs(Th - ((T0 + T1)/2)) < 0.0001:
                            break
                        else:
                            Th = (T0 + T1)/2 
                            
            elif method == "otsu":
                for ii in trange(len(gene_data)):
                    fig1 = gene_data[ii, :]
                    r1 = raw_gene_data[ii, :] 
                    thres_list = np.arange(0.01,0.995,0.025)
                    temp_std = np.zeros(thres_list.shape)
                    for iii in range(len(thres_list)):
                        temp_thres = thres_list[iii]
                        q1 = fig1 > temp_thres
                        b1 = fig1 <= temp_thres
                        qv = r1[q1]
                        bv = r1[b1]
                        if len(qv) >= len(r1) * 0.15:
                            temp_std[iii] = (len(qv) * np.std(qv) + len(bv) * np.std(bv)) / gene_data.shape[1]
                        else:
                            temp_std[iii] = 1e4
                    Th = thres_list[temp_std == np.min(temp_std)]
                    
            # a1 = fig1 >= Th
            # I, _ = make_image(a1, locates)
            output[i, :] = fig1 >= Th
    #--------------------------------------------------------------------------        
    else:               
        output = np.zeros((gene_data.shape[0],adata.uns['shape'][0]*adata.uns['shape'][1]))  
        for i in trange(len(gene_data)):
            fig1 = gene_data[i, :]

            if method == "iterative":    
                    # m, n = Im.shape       
                    zd = float(np.nanmax(fig1))
                    zx = float(np.nanmin(fig1))
                    Th = float((zd+zx))/2
                    while True:                       
                        S1 = np.sum(fig1[fig1>=Th])
                        n1 = len(fig1[fig1>=Th])
                        S0 = np.sum(fig1[fig1<Th])
                        n0 = len(fig1[fig1<Th])
                        T0 = S0/n0; T1 = S1/n1
                                        
                        if abs(Th - ((T0 + T1)/2)) < 0.0001:
                            break
                        else:
                            Th = (T0 + T1)/2 
                            
            elif method == "otsu":
                for ii in trange(len(gene_data)):
                    fig1 = gene_data[ii, :]
                    
                    img = fig1.reshape(adata.uns['shape'])
                    Th2, a_img = cv2.threshold(img.astype(np.uint8), 0, 255, cv2.THRESH_OTSU)

                    
            a1 = fig1 >= Th
            # I, _ = make_image(a1, locates)
            output[i, :] = fig1 >= Th
    #--------------------------------------------------------------------------        
    adata.uns['binary_image'] = output+0
    return adata


def get_sub(adata, connect_kernel_size = 5, neighbors = 8, platform="visium",del_rate = 0.01): 
    gene_data = adata.uns['binary_image']
    locates = adata.uns['locates']
        
    print('\nSpliting subregions for each gene:')
    #--------------------------------------------------------------------------
    if platform=="visium":
        image_idx_1d = adata.obs['image_idx_1d']
        output = np.zeros(gene_data.shape)
        del_index = np.ones(gene_data.shape[0])
        for i in trange(len(gene_data)):
            temp_data = gene_data[i, :]
            temp_i, _ = make_image(temp_data, locates)      
            kernel = np.ones((connect_kernel_size,connect_kernel_size), np.uint8)
            temp_i = cv2.morphologyEx(temp_i, cv2.MORPH_CLOSE, kernel) # close
            # plt.imshow(temp_i)
            region_label = label(temp_i, neighbors) - 1
            T = np.zeros(region_label.shape)
            classes = np.max(np.unique(region_label)) + 1      
            len_list = np.zeros(classes)     
            for j in range(classes):
                len_list[j] = len(region_label[region_label == j])
            cond = len_list >= gene_data.shape[1] * 0.01        
            if len(np.where(cond[1:] == True)[0]) == 0:
                # tqdm.write("\nThe {}th gene will be deleted".format(i+1))
                del_index[i] = 0
                # continue
            indexes = np.where(cond == True)[0]       
            for j in range(len(indexes)):
                tar_num = indexes[j]
                tar_locs = region_label == tar_num
                T[tar_locs] = j
            targe_image = T * (temp_i > 0)
            classes_n = np.max(np.unique(targe_image)).astype(int) + 1       
            len_list_n = np.zeros(classes_n)        
            for j in range(classes_n):
                len_list_n[j] = len(targe_image[targe_image == j])            
            if len(len_list_n) > 1:                        
                if np.max(len_list_n[1:]) < gene_data.shape[1] * del_rate:
                    # tqdm.write("\nThe {}th gene will be deleted".format(i+1))
                    del_index[i] = 0
                    # continue
            else:
                del_index[i] = 0
                # continue
            output[i, :] = gene_img_flatten(targe_image, image_idx_1d)      
            # adata.uns['subregions'] = output
            # adata.uns['del_index'] = del_index.astype(int)
    #--------------------------------------------------------------------------
    else:
        output = np.zeros((gene_data.shape[0], adata.uns['shape'][0]*adata.uns['shape'][1]))
        del_index = np.ones(gene_data.shape[0])
        for i in trange(len(gene_data)):
            temp_data = gene_data[i, :]
            temp_i = temp_data.reshape(adata.uns['shape'])     
            kernel = np.ones((connect_kernel_size,connect_kernel_size), np.uint8)
            temp_i = cv2.morphologyEx(temp_i, cv2.MORPH_CLOSE, kernel)
            # plt.imshow(temp_i)
            region_label = label(temp_i, neighbors) - 1
            #plt.imshow(region_label)
            T = np.zeros(region_label.shape)
            classes = np.max(region_label) + 1      
            len_list = np.zeros(classes)
            for j in range(classes):
                len_list[j] = len(region_label[region_label == j])
            cond = len_list >= gene_data.shape[1] * 0.002        
            if len(np.where(cond[1:] == True)[0]) == 0:
                # tqdm.write("\nThe {}th gene will be deleted".format(i+1))
                del_index[i] = 0
                # continue
            indexes = np.where(cond == True)[0]       
            for j in range(len(indexes)):
                tar_num = indexes[j]
                tar_locs = region_label == tar_num
                T[tar_locs] = j
            targe_image = T * (temp_i > 0)
            # plt.imshow(targe_image )
            classes_n = np.max(np.unique(targe_image)).astype(int) + 1       
            len_list_n = np.zeros(classes_n)        
            for j in range(classes_n):
                len_list_n[j] = len(targe_image[targe_image == j])            
            if len(len_list_n) > 1:                        
                if np.max(len_list_n[1:]) < gene_data.shape[1] * del_rate:
                    # tqdm.write("\nThe {}th gene will be deleted".format(i+1))
                    del_index[i] = 0
                    # continue
            else:
                del_index[i] = 0
                # continue
            output[i, :] = targe_image.flatten()
    #--------------------------------------------------------------------------
    adata.uns['subregions'] = output
    adata.uns['del_index'] = del_index.astype(int)  
    return adata


def cal_PI(adata,platform="visium"):
    data = adata.uns['nor_counts']
    subregions = adata.uns['subregions']
    del_idx = adata.uns['del_index']
    
    print('\nComputing PROST Index for each gene:')
    #--------------------------------------------------------------------------
    if platform=="visium":
        SEP = np.zeros(len(data))
        SIG = np.zeros(len(data))
        region_number = np.zeros(len(data))
        for i in trange(len(data)): 
            temp_raw = data[i, :]
            temp_label = subregions[i, :]
            back_value = temp_raw[temp_label == 0]
            back_value = back_value[back_value > 0]
            
            if back_value.size == 0:
                back_value = 0  
            class_mean = np.zeros(max(np.unique(temp_label)).astype(int) + 1)
            class_var = np.zeros(max(np.unique(temp_label)).astype(int) + 1)
            class_std = np.zeros(max(np.unique(temp_label)).astype(int) + 1)
            class_len = np.zeros(max(np.unique(temp_label)).astype(int) + 1)
         
            for ii in range(max(np.unique(temp_label)).astype(int) + 1):
                Temp = temp_raw[temp_label == ii]
                if Temp.size == 0:
                    class_value = 0
                else:
                    class_value = Temp
                class_mean[ii] = np.mean(class_value)
                class_var[ii] = np.var(class_value)
                class_std[ii] = np.std(class_value)
                if isinstance(class_value, int):
                    if class_value == 0:
                        class_len[ii] = 0
                    else:              
                        class_len[ii] = len(class_value) - 1               
                else:
                    if class_value.size == 0:
                        class_len[ii] = 0
                    else:              
                        class_len[ii] = len(class_value) - 1
                        
            target_class = np.where(class_mean > 0)[0]
            class_mean = class_mean[target_class]
            class_std = class_std[target_class]
            class_var = class_var[target_class]
            class_len = class_len[target_class]
            
            # Calculate Separability and Significance
            SEP[i] = 1 - sum((class_len * class_var)) / ((len(temp_raw)-1) * np.var(temp_raw))
            SIG[i] = (np.mean(class_mean) - np.mean(back_value)) / sum(class_std / class_mean) 
            region_number[i] = len(class_len)
            del class_mean, class_var, class_len, class_std
            
        # Pattern Index    
        PI = minmax_normalize(SEP) * minmax_normalize(SIG)
        PI = PI * del_idx
        adata.var["SEP"] = SEP
        adata.var["SIG"] = SIG
        adata.var["PI"] = PI
        adata.var["n_regions"] = region_number
    #--------------------------------------------------------------------------
    else:
        locates = adata.uns['locates']
        SEP = np.zeros(len(data))
        SIG = np.zeros(len(data))
        for i in trange(len(data)): 
            temp_raw = data[i, :]
            temp_img,_ = make_image(temp_raw, locates, platform)
            temp_raw = temp_img.flatten()
            temp_label = subregions[i, :]
            back_value = temp_raw[temp_label == 0]
            back_value = back_value[back_value > 0]
            
            if back_value.size == 0:
                back_value = 0  
            class_mean = np.zeros(max(np.unique(temp_label)).astype(int) + 1)
            class_var = np.zeros(max(np.unique(temp_label)).astype(int) + 1)
            class_std = np.zeros(max(np.unique(temp_label)).astype(int) + 1)
            class_len = np.zeros(max(np.unique(temp_label)).astype(int) + 1)
         
            for ii in range(max(np.unique(temp_label)).astype(int) + 1):
                Temp = temp_raw[temp_label == ii]
                if Temp.size == 0:
                    class_value = 0
                else:
                    class_value = Temp
                class_mean[ii] = np.nanmean(class_value)
                class_var[ii] = np.nanvar(class_value)
                class_std[ii] = np.nanstd(class_value)
                if isinstance(class_value, int):
                    if class_value == 0:
                        class_len[ii] = 0
                    else:              
                        class_len[ii] = len(class_value) - 1               
                else:
                    if class_value.size == 0:
                        class_len[ii] = 0
                    else:              
                        class_len[ii] = len(class_value) - 1
                        
            target_class = np.where(class_mean > 0)[0]
            class_mean = class_mean[target_class]
            class_std = class_std[target_class]
            class_var = class_var[target_class]
            class_len = class_len[target_class]
            
            # Calculate Separability and Significance
            SEP[i] = 1 - sum((class_len * class_var)) / ((len(temp_raw)-1) * np.nanvar(temp_raw))
            SIG[i] = (np.mean(class_mean) - np.mean(back_value)) / sum(class_std / class_mean)      
            del class_mean, class_var, class_len, class_std
            
        # Pattern Index    
        PI = minmax_normalize(SEP) * minmax_normalize(SIG)
        PI = PI * del_idx
        adata.var["SEP"] = SEP
        adata.var["SIG"] = SIG
        adata.var["PI"] = PI
        
        adata.uns['shape'] = []
    #--------------------------------------------------------------------------
    adata.uns['nor_counts'] = []
    adata.uns['binary_image'] = []
    adata.uns['subregions'] = []
    adata.uns['del_index'] = []
    
    return adata
