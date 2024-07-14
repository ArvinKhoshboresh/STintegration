import anndata
import sys
import scanpy as sc
import numpy as np

np.set_printoptions(edgeitems=30)

# Load the h5ad file
h5ad_path = sys.argv[1]
adata = sc.read_h5ad(h5ad_path)

print(adata)

print(adata.obs_names)

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
