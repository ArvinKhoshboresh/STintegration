import anndata
import sys
import scanpy as sc

# Load the h5ad file
h5ad_path = sys.argv[1]
adata = sc.read_h5ad(h5ad_path)

print(adata)

# Extract spatial coordinates
coords1 = adata.obsm['X_spatial']

# Print the gene matrix
print("Coords (X_spatial):")
print(adata.X)

# Create matrix of gene expression x cell
expression_matrix1 = adata.X

# Print the gene matrix
print("Gene Matrix (adata.X):")
print(adata.X)

# Print cluster annotations if they exist
if 'clusters' in adata.obs:
    print("\nCluster Annotations (adata.obs['clusters']):")
    print(adata.obs['clusters'])
else:
    print("\nCluster Annotations not found in adata.obs")

# Print spatial coordinates if they exist
if 'spatial' in adata.obsm:
    print("\nSpatial Coordinates (adata.obsm['spatial']):")
    print(adata.obsm['spatial'])
else:
    print("\nSpatial Coordinates not found in adata.obsm")

# For spatial coordinates with xyz
if 'xyz' in adata.obsm:
    print("\nXYZ Spatial Coordinates (adata.obsm['xyz']):")
    print(adata.obsm['xyz'])
else:
    print("\nXYZ Spatial Coordinates not found in adata.obsm")
