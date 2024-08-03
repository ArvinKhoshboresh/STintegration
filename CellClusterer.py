import scanpy as sc
import anndata as ad
import numpy as np
from sklearn.cluster import KMeans


def combine_and_cluster(adatas, n_clusters=10):
    # Combine multiple AnnData objects
    combined_data = ad.concat(adatas, axis=0)

    # Calculate Euclidean distances and perform clustering
    X = combined_data.X.toarray() if hasattr(combined_data.X, 'toarray') else combined_data.X
    kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(X)

    start_idx = 0
    for adata in adatas:
        end_idx = start_idx + adata.n_obs
        adata.obs['Kclusters'] = kmeans.labels_[start_idx:end_idx].astype(str)
        start_idx = end_idx

    return adatas


files = [
    "/home/arvin/barseq/controls/Control_1.mat.h5ad",
    "/home/arvin/barseq/controls/Control_2.mat.h5ad",
    "/home/arvin/barseq/controls/Control_3.mat.h5ad",
    "/home/arvin/barseq/controls/Control_4.mat.h5ad",
    "/home/arvin/barseq/enucleated/Enucleated_1.mat.h5ad",
    "/home/arvin/barseq/enucleated/Enucleated_2.mat.h5ad",
    "/home/arvin/barseq/enucleated/Enucleated_3.mat.h5ad",
    "/home/arvin/barseq/enucleated/Enucleated_4.mat.h5ad"]
adata_struct = [sc.read(file) for file in files]

# Combine and cluster
adatas = combine_and_cluster(adata_struct, n_clusters=200)

for idx, adata in enumerate(adatas):
    adata.write_h5ad(filename=f"{files[idx].split('.')[0]}-DeNovoClustered.h5ad")
