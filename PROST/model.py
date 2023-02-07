import pandas as pd
import numpy as np
import scanpy as sc
import torch
import torch.nn as nn
# import torch.nn.functional as F
from torch.nn.parameter import Parameter
import torch.optim as optim
# device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
# from sklearn.mixture import GaussianMixture
# from sklearn.mixture import BayesianGaussianMixture
from scipy.sparse import issparse
from tqdm import tqdm
from . layers import *
from . utils import *



class PROST_NN(nn.Module):
    def __init__(self, nfeat, nhid, dropout, leaky_alpha, beta = 0.5):
        super(PROST_NN, self).__init__()
        self.gal = GraphAttentionLayer(nfeat, nhid, dropout, leaky_alpha)
        self.nhid=nhid
        # self.nclass=nclass
        self.beta=beta

    def forward(self, x, adj):
        x = self.gal(x, adj)
        q = 1.0 / ((1.0 + torch.sum((x.unsqueeze(1) - self.mu)**2, dim=2) / self.beta) + 1e-8)
        q = q**(self.beta+1.0)/2.0
        q = q / torch.sum(q, dim=1, keepdim=True)
        return x, q

    def loss_function(self, p, q):
        loss = torch.mean(torch.sum(p*torch.log(p/(q+1e-6)), dim=1))
        return loss

    def target_distribution(self, q):
        p = q**2 / torch.sum(q, dim=0)
        p = p / torch.sum(p, dim=1, keepdim=True)
        return p

    def fit(self, 
            X, 
            adj,  
            lr=0.001, 
            max_epochs=5000, 
            update_interval=3, 
            weight_decay=5e-4,
            init="louvain",
            n_neighbors=10,
            res=0.4,
            n_clusters=7,
            tol=1e-3,
            seed = 818):

        optimizer = optim.Adam(self.parameters(),lr=lr, weight_decay=weight_decay)

        features = self.gal(torch.FloatTensor(X),torch.FloatTensor(adj))
     
        #----------------------------------------------------------------           
        if init=="kmeans":
            print("\nInitializing cluster centers with kmeans, n_clusters known")
            self.n_clusters=n_clusters
            kmeans = KMeans(self.n_clusters, n_init=20)
            y_pred = kmeans.fit_predict(features.detach().numpy())
            
        elif init=="mclust":
            print("\nInitializing cluster centers with mclust, n_clusters known")
            data = features.detach().numpy()
            self.n_clusters = n_clusters
            self.seed = seed
            y_pred = mclust_R(data, num_cluster = self.n_clusters, random_seed = self.seed)
            # y_pred = y_pred.dropna()
            y_pred = y_pred.astype(int)

        elif init=="louvain":
            print("\nInitializing cluster centers with louvain, resolution = ", res)
            adata = sc.AnnData(features.detach().numpy())
            sc.pp.neighbors(adata, n_neighbors=n_neighbors)
            sc.tl.louvain(adata, resolution=res)
            y_pred = adata.obs['louvain'].astype(int).to_numpy()
            self.n_clusters = len(np.unique(y_pred))
            
        elif init=="leiden":
            print("\nInitializing cluster centers with leiden, resolution = ", res)
            adata=sc.AnnData(features.detach().numpy())
            sc.pp.neighbors(adata, n_neighbors=n_neighbors)
            sc.tl.leiden(adata, resolution=res)
            y_pred = adata.obs['leiden'].astype(int).to_numpy()
            self.n_clusters = len(np.unique(y_pred))
        #----------------------------------------------------------------
        y_pred_last = y_pred
        self.mu = Parameter(torch.Tensor(self.n_clusters, self.nhid))
        X = torch.FloatTensor(X)
        adj = torch.FloatTensor(adj)


        features = pd.DataFrame(features.detach().numpy()).reset_index(drop = True)
        Group = pd.Series(y_pred, index=np.arange(0,features.shape[0]), name="Group")
        Mergefeature = pd.concat([features,Group], axis=1)
        cluster_centers = np.asarray(Mergefeature.groupby("Group").mean())       
        self.mu.data.copy_(torch.Tensor(cluster_centers))
        
        #----------------------------------------------------------------
        self.train()
        with tqdm(total=max_epochs) as t:
            for epoch in range(max_epochs):
                
                t.set_description('Epoch')
                if epoch%update_interval == 0:
                    _, q = self.forward(X,adj)
                    p = self.target_distribution(q).data
                    t.update(update_interval)
                optimizer.zero_grad()
                z,q = self(X, adj)
                loss = self.loss_function(p, q)
                loss.backward()
                optimizer.step()
    
                # if epoch%10==0:
                #     print("Epoch:", epoch,"\tLoss:",loss.data.numpy())
                t.set_postfix(loss = loss.data.numpy())
                
                #Check stop criterion
                y_pred = torch.argmax(q, dim=1).data.cpu().numpy()
                delta_label = np.sum(y_pred != y_pred_last).astype(np.float32) / X.shape[0]
                y_pred_last = y_pred
                if epoch>0 and (epoch-1)%update_interval == 0 and delta_label < tol:
                    print('delta_label ', delta_label, '< tol ', tol)
                    print("Reach tolerance threshold. Stopping training.")
                    print("Total epoch:", epoch)
                    break

    def predict(self, X, adj):
        z,q = self(torch.FloatTensor(X),torch.FloatTensor(adj))
        return z, q
    
    
class PROST_cluster(object):
    def __init__(self):
        super(PROST_cluster, self).__init__()

    def train(self,
                embed,
                adj, 
                n_pcs=15, 
                lr=0.005,
                dropout=0.2,
                leaky_alpha=0.2,
                max_epochs=2000,
                weight_decay=5e-4,
                init="louvain", #louvain or leiden
                n_neighbors=10, #for louvain or leiden
                n_clusters=7, #for mclust or kmeans
                res=0.4, #for louvain or leiden
                seed = 818, #for mclust
                tol=1e-3):
        self.n_pcs=n_pcs
        self.res=res
        self.lr=lr
        self.dropout=dropout
        self.leaky_alpha=leaky_alpha
        self.max_epochs=max_epochs
        self.weight_decay=weight_decay
        self.init=init
        self.n_neighbors=n_neighbors
        self.n_clusters=n_clusters
        self.res=res
        self.tol=tol
        self.seed=seed
        # assert adata.shape[0]==adj.shape[0]==adj.shape[1]
        
        # pca = PCA(n_components=self.n_pcs)
        # print("\nRunning PCA ...")
        if issparse(embed):
            embed=embed.A
        # else:
        #     embed=embed
         
        self.model = PROST_NN(embed.shape[1], embed.shape[1], self.dropout, self.leaky_alpha)

        self.model.fit(embed,
                       adj,
                       lr=self.lr,
                       max_epochs=self.max_epochs,
                       weight_decay=self.weight_decay,
                       init=self.init,
                       n_neighbors=self.n_neighbors,
                       n_clusters=self.n_clusters,
                       res=self.res, 
                       tol=self.tol,
                       seed=self.seed)
        
        self.embed = embed
        self.adj = adj

    def predict(self):
        z,q = self.model.predict(self.embed, self.adj)
        y_pred = torch.argmax(q, dim=1).data.cpu().numpy()
        return y_pred, q.detach().numpy(), z.detach().numpy()