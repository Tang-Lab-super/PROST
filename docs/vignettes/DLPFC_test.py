import pandas as pd
import numpy as np
import scanpy as sc
import os
import warnings
warnings.filterwarnings("ignore")
import matplotlib as mpl
import matplotlib.pyplot as plt

import PROST


# the location of R (used for the mclust clustering)
ENVpath = "your path of PROST_ENV"            # refer to 'How to use PROST' section
os.environ['R_HOME'] = f'{ENVpath}/lib/R'
os.environ['R_USER'] = f'{ENVpath}/lib/python3.7/site-packages/rpy2'

#%% Set seed
SEED = 818
PROST.setup_seed(SEED)

#%% Read in data
section_num = 151672

# Set directory (If you want to use additional data, please change the file path)
rootdir = 'datasets/DLPFC'

input_dir = os.path.join(f'{rootdir}', str(section_num))
spatial_dir = os.path.join(f'{rootdir}', str(section_num),'spatial')
output_dir = os.path.join(f'{rootdir}', str(section_num), 'results')
if not os.path.isdir(output_dir):
    os.makedirs(output_dir)
    
#
adata = sc.read_visium(path=input_dir, count_file='{}_filtered_feature_bc_matrix.h5'.format(section_num))
adata.var_names_make_unique()

#%% Calculate PI
adata = PROST.prepare_for_PI(adata, platform="visium")
adata = PROST.cal_prost_index(adata, platform="visium")

# Save PI result
adata.write_h5ad(output_dir+"/PI_result.h5")

# Draw SVGs detected by PI
PROST.plot_gene(adata, 
                platform="visium",
                size = 2, 
                sorted_by = "PI", 
                top_n = 50, 
                save_path = output_dir)

# Calculate Moran'I and Geary'C for SVGs dected by PI
PROST.cal_moran_I_and_geary_C_for_PI_SVGs(adata, PI_top_n=50, save_path = output_dir)

#%% Clustering
# Set the number of clusters
n_clusters = 5

# 1.Read PI result
adata = sc.read(output_dir+"/PI_result.h5")

# 2.Expression data preprocessing
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata = PROST.feature_selection(adata, save_path = output_dir, by = "prost", n_top_genes = 3000)

# 3.Run PROST clustering
PROST.run_prost_clust(adata,
                      platform="visium",
                      key_added = "PROST",
                      init="mclust",                         
                      n_clusters = n_clusters,                        
                      gnnlayers = 2,                                                                       
                      laplacin_filter = True,                        
                      lr = 0.1,                         
                      SEED=SEED,                          
                      max_epochs = 500,                        
                      tol = 5e-3,                        
                      post_processing = True,                        
                      pp_run_times = 3)

adata.write_h5ad(output_dir+"/PNN_result.h5")   

# 4.Save result
clustering = adata.obs["clustering"]
clustering.to_csv(output_dir+"/clusters.csv",header = False)
pp_clustering = adata.obs["pp_clustering"] 
pp_clustering.to_csv(output_dir+"/pp_clusters.csv",header = False)

#%% Plot preparation
# Read data
adata = sc.read(output_dir+'/PNN_result.h5')

# Read annotation
labels_true = pd.read_csv(input_dir+'/cluster_labels.csv')
labels_true.index = labels_true["key"].str[7:]
adata.obs["annotation"] = labels_true["ground_truth"]
adata.obs["annotation"] = adata.obs["annotation"].astype('category').astype('str')

# Set colors
plot_color = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2',
              "#F56867","#FEB915","#C798EE","#59BE86","#7495D3","#D1D1D1","#6D1A9C",
              "#15821E","#3A84E6","#997273","#DB4C6C","#554236",
              "#AF5F3C","#93796C","#F9BD3F","#DAB370","#268785"]
cmp = mpl.colors.ListedColormap(plot_color)

#%% Plot cluster result
# read ground_truth
labels_true = pd.read_csv(input_dir+'/cluster_labels.csv')
labels_true = labels_true.fillna("NA")
adata.obs["annotation"] = np.array(labels_true['ground_truth'])

# plot
plt.rcParams["figure.figsize"] = (4,4)
ax = sc.pl.spatial(adata, 
                  img_key = "hires", 
                  color = ["annotation","clustering","pp_clustering"],
                  title = ["Manual annotation",'clustering','post_processed clustering'],
                  color_map = cmp, 
                  na_in_legend = False,
                  ncols = 3,
                  size = 1.1,
                  # wspace = 0.1,
                  show = False)
plt.savefig(output_dir+"/cluster_result.png", dpi=600)
plt.close()

#%% Calculate ARI and NMI
ARI, NMI, silhouette_score = PROST.cal_metrics_for_DLPFC(adata.obs["pp_clustering"], labels_true_path = input_dir+'/cluster_labels.csv')
#save result
single_data = {'ARI':[ARI], 'NMI':[NMI], 'silhouette_score':[silhouette_score]}
metrics = pd.DataFrame(single_data)
metrics.to_csv(output_dir+'/clustering_metrics.csv', index = False)

#%% Plot PAGA
plt.rcParams["figure.figsize"] = (4,3)
used_adata = adata[adata.obs["annotation"]!='nan']
sc.pp.neighbors(used_adata, use_rep="PROST")
sc.tl.umap(used_adata)
sc.tl.paga(used_adata,groups='annotation')
ax = sc.pl.paga_compare(used_adata, color="annotation",frameon=False, show = False, random_state = SEED,fontoutline = 2)
plt.savefig(output_dir+"/PAGA_compare.png", dpi=600)
plt.close()
