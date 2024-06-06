import scanpy as sc
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import sys

# Load the MERFISH data from an h5ad file
h5ad_path = sys.argv[1]
adataFull = sc.read_h5ad(h5ad_path)
print('data loaded')
print(adataFull)

#################### UMAPPING ###############

print('started Umapping')

# Inspect the first few rows and columns of the data matrix .X
print("First few entries in the data matrix (genes in rows):")
print(adataFull.X[:5, :5].toarray() if hasattr(adataFull.X, "toarray") else adataFull.X[:5, :5])

# Inspect the first few rows of the cell metadata
print("Cell metadata (.obs):")
print(adataFull.obs.head())

# Inspect the first few rows of the gene metadata
print("Gene metadata (.var):")
print(adataFull.var.head())

# # Cut data into pieces for faster prototyping
# num_cells = adataFull.shape[0]
# indices = np.random.permutation(num_cells)
# split = indices[:num_cells // 15]
# adata = adataFull[split].copy()
# print('data split')

adata = adataFull.copy()

# Preprocess the data: Normalize and log-transform if needed
sc.pp.normalize_total(adata, target_sum=1e6, exclude_highly_expressed=True, inplace=True)
#sc.pp.log1p(adata)
print('preprocessing done')

# Run PCA to reduce dimensionality
sc.pp.pca(adata, n_comps=20)
print('PCA done')
print(adata)

# Compute neighbors and UMAP on PCA-reduced data
sc.pp.neighbors(adata, use_rep='X_pca')
print('Neighbours calced')
print(adata)
sc.tl.umap(adata)
print('umap calced')
print(adata)

# adata.write("C:/Users/arkho/OneDrive/Desktop/pilot_processed_adata.h5ad")

# # Plot UMAP
# with plt.rc_context({"figure.figsize": (12, 12), "figure.dpi": (500)}):
#     sc.pl.umap(adata, save='umap_plot.png')
# print('umap plotted')

##################### SPATIAL MAPPING ##################

print('started spatial mapping')

# Load spatial coordinates from the CSV file
coords_path = sys.argv[2]
spatial_coords = pd.read_csv(coords_path)

# Ensure the indices match between the spatial coordinates and adata
spatial_coords.set_index('cell_label', inplace=True)
print(spatial_coords)


# Find the common cell labels
common_cells = adata.obs_names.intersection(spatial_coords.index)

# Subset the AnnData object to keep only common cells
adata = adata[common_cells].copy()

# Subset the spatial coordinates to keep only common cells
spatial_coords = spatial_coords.loc[common_cells]

# Check the dimensions after filtering
print(adata)
print(spatial_coords.shape)


# Check that the cell labels match
if not all(spatial_coords.index.isin(adata.obs_names)):
    raise ValueError("Some cell labels in the spatial coordinates do not match those in the AnnData object.")
if not all(adata.obs_names.isin(spatial_coords.index)):
    raise ValueError("Some cell labels in the AnnData object do not have corresponding spatial coordinates.")

# Reorder spatial_coords to match the order of adata.obs_names
spatial_coords = spatial_coords.loc[adata.obs_names]

# Add spatial coordinates to the AnnData object
adata.obsm['X_spatial'] = spatial_coords[['x', 'y', 'z']].values

print(adata.obsm['X_spatial'])
print('calculation end')

# # Plot the spatial distribution
# with plt.rc_context({"figure.figsize": (12, 12), "figure.dpi": (500)}):
#     sc.pl.spatial(adata, spot_size=0.01, save='spatial_plot.png')
# print('spatial plotted')

adata.write("pilot_adata_UMAPed_Spatialed.h5ad")

print('script ended')
