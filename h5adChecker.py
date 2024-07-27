import anndata
import sys
import scanpy as sc
import numpy as np
from matplotlib import pyplot as plt

np.set_printoptions(edgeitems=30)

# Load the h5ad file
h5ad_path = sys.argv[1]
h5ad_path = sys.argv[2]
h5ad_path = sys.argv[3]
h5ad_path = sys.argv[4]
adata1 = sc.read_h5ad(h5ad_path)
adata2 = sc.read_h5ad(h5ad_path)
adata3 = sc.read_h5ad(h5ad_path)
adata4 = sc.read_h5ad(h5ad_path)

print(adata1)

print(adata1.obs["subclass"].nunique())


def extract_coords(adata):
    """
    Extracts the CCF_AP_axis, CCF_ML_axis, and CCF_DV_axis columns from the obs DataFrame of an anndata object
    and combines them into a single matrix.

    Parameters:
    adata (anndata.AnnData): The anndata object.

    Returns:
    np.ndarray: A matrix with X, Y, and Z coordinates in each row.
    """
    x = adata.obs['CCF_AP_axis'].values
    y = adata.obs['CCF_ML_axis'].values
    z = adata.obs['CCF_DV_axis'].values
    return np.vstack((x, y, z)).T


def flip_coordinates(data):
    # Assuming data is a numpy array of shape (n, 3) with columns [x, y, z]
    flipped_data = data.copy()
    flipped_data[:, 0] = flipped_data[:, 0] * -1  # Flip x-coordinates
    return flipped_data

# coords1 = extract_coords(adata1)
# coords1 = flip_coordinates(coords1)
#
# adata1.obs['CCF_AP_axis'].values = coords1[0]
# adata1.obs['CCF_ML_axis'].values = coords1[1]
# adata1.obs['CCF_DV_axis'].values = coords1[2]

# Cut data into pieces for faster prototyping
cut_data = True
if cut_data:
    cut_data_factor = 80
    num_cells = adata1.shape[0]
    indices = np.random.permutation(num_cells)
    split = indices[:num_cells // cut_data_factor]
    adata1 = adata1[split].copy()

    num_cells = adata2.shape[0]
    indices = np.random.permutation(num_cells)
    split = indices[:num_cells // cut_data_factor]
    adata2 = adata2[split].copy()

    num_cells = adata3.shape[0]
    indices = np.random.permutation(num_cells)
    split = indices[:num_cells // cut_data_factor]
    adata3 = adata3[split].copy()

    num_cells = adata4.shape[0]
    indices = np.random.permutation(num_cells)
    split = indices[:num_cells // cut_data_factor]
    adata4 = adata4[split].copy()

    print('WARNING: Data split')

coords1 = extract_coords(adata1)
coords1 = flip_coordinates(coords1)
coords2 = extract_coords(adata2)
coords3 = extract_coords(adata3)
coords4 = extract_coords(adata4)

# Plot this region in Coronal view
fig, neighours_fig = plt.subplots(figsize=(15, 10), dpi=800)
# Plot all cells from both datasets
neighours_fig.scatter(coords1[:, 2], coords1[:, 1], c='blue', label='adata1', alpha=0.5, s=0.2, linewidths=0.3)
neighours_fig.scatter(coords2[:, 2], coords2[:, 1], c='red', label='adata2', alpha=0.5, s=0.2, linewidths=0.3)
neighours_fig.scatter(coords3[:, 2], coords3[:, 1], c='green', label='adata3', alpha=0.5, s=0.2, linewidths=0.3)
neighours_fig.scatter(coords4[:, 2], coords4[:, 1], c='orange', label='adata4', alpha=0.5, s=0.2, linewidths=0.3)
# for cell_Idx, cell_coord in enumerate(coords1):
#     neighours_fig.annotate(adata1.obs['slice'][cell_Idx], (cell_coord[2], cell_coord[1]), fontsize=1)
# for cell_Idx, cell_coord in enumerate(coords2):
#     neighours_fig.annotate(adata2.obs['slice'][cell_Idx], (cell_coord[2], cell_coord[1]), fontsize=1)
neighours_fig.set_xlabel('Z Coordinate')
neighours_fig.set_ylabel('Y Coordinate')
plt.savefig('graph.png', bbox_inches='tight')

# # Extract the specified columns from obs
# columns_to_extract = ['CCFano', 'CCFname', 'clustid', 'clustname', 'subclass']
#
# # Ensure the columns exist in the obs DataFrame
# for col in columns_to_extract:
#     if col not in adata.obs.columns:
#         raise KeyError(f"Column '{col}' not found in the AnnData object's obs attribute.")
#
# # Print the values for each row
# for index, row in adata.obs[columns_to_extract].iterrows():
#     print(f"Row {index}: {row.to_dict()}")
#
# # Create a list of all unique CCFano entries
# unique_ccfano = adata.obs['CCFname'].unique().tolist()
#
# # Print the unique CCFano entries
# print("Unique CCFano entries:")
# for entry in unique_ccfano:
#     print(entry)

# # Print the spatial data
# print(f"Gene Matrix Type : {type(adata.obs.CCF_AP_axis)}")
# print("ap axis (adata.obs.CCF_AP_axis):")
# print(adata.obs.CCF_AP_axis)
#
# # Print the spatial data
# print(f"Gene Matrix Type : {type(adata.obs.CCF_ML_axis)}")
# print("ap axis (adata.obs.CCF_ML_axis):")
# print(adata.obs.CCF_ML_axis)
#
# # Print the spatial data
# print(f"Gene Matrix Type : {type(adata.obs.CCF_DV_axis)}")
# print("ap axis (adata.obs.CCF_DV_axis):")
# print(adata.obs.CCF_DV_axis)
#
# # Extract the columns from the 'obs' DataFrame
# x = adata.obs['CCF_AP_axis'].values
# y = adata.obs['CCF_ML_axis'].values
# z = adata.obs['CCF_DV_axis'].values
#
# print(f"X- min:{min(x)} max:{max(x)} ")
# print(f"Y- min:{min(y)} max:{max(y)} ")
# print(f"Z- min:{min(z)} max:{max(z)} ")

# # Combine these columns into a single matrix
# xyz_matrix = np.vstack((x, y, z)).T
#
# print(xyz_matrix)

# adata.obs.CCF_AP_axis.to_csv("series", sep=',', encoding='utf-8')

# # Print the gene matrix
# print(f"Gene Matrix Type : {type(adata.X)}")
# print("Gene Matrix (adata.X):")
# print(adata.X.todense())
#
# # Print the gene list
# print(f"Gene Matrix Type : {type(adata.var)}")
# print("Gene Matrix (adata.var):")
# print(adata.var)
