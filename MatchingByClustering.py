import numpy as np
from sklearn.cluster import KMeans
import sys
import scanpy as sc
from scipy.optimize import linear_sum_assignment

adata1_path = sys.argv[1]
adata2_path = sys.argv[2]
adata3_path = sys.argv[3]
adata4_path = sys.argv[4]

adata_struct = []
adata_struct.append(sc.read_h5ad(adata1_path))
adata_struct.append(sc.read_h5ad(adata2_path))
adata_struct.append(sc.read_h5ad(adata3_path))
adata_struct.append(sc.read_h5ad(adata4_path))

num_neighbours = 150

cut_data = True
np.random.seed(42)  # For reproducibility
if cut_data:
    cut_data_factor = 100000
    for idx in range(0, len(adata_struct)):
        num_cells = adata_struct[idx].shape[0]
        indices = np.random.permutation(num_cells)
        split = indices[:num_cells // cut_data_factor]
        adata_struct[idx] = adata_struct[idx][split].copy()

    print('WARNING: Data split')


def extract_coords(adata):
    x = adata.obs['CCF_AP_axis'].values
    y = adata.obs['CCF_ML_axis'].values
    z = adata.obs['CCF_DV_axis'].values
    return np.vstack((x, y, z)).T


# Extract spatial coordinates
coords1 = extract_coords(adata_struct[0])
coords2 = extract_coords(adata_struct[1])
coords3 = extract_coords(adata_struct[2])
coords4 = extract_coords(adata_struct[3])

datasets = [coords1, coords2, coords3, coords4]
# Determine the smallest dataset size
min_size = min(len(data) for data in datasets)

# Concatenate datasets for clustering
all_data = np.concatenate(datasets, axis=0)
print(f"All data: {all_data}\n")

# Fit k-means to form clusters
kmeans = KMeans(n_clusters=min_size).fit(all_data)
cluster_centers = kmeans.cluster_centers_

# Create a cost matrix for the assignment problem
cost_matrices = [np.zeros((len(dataset), min_size)) for dataset in datasets]

tracker = 0
for i, data in enumerate(datasets):
    start_idx = tracker
    end_idx = start_idx + len(data)
    # Compute pairwise distances
    distances = np.linalg.norm(data[:, np.newaxis] - cluster_centers, axis=2)
    cost_matrices[i] = distances

print(cost_matrices)

# Solve the assignment problem
row_ind, col_ind = linear_sum_assignment(cost_matrices)

# Form the final clusters
final_clusters = [[] for _ in range(min_size)]
for i in range(len(row_ind)):
    dataset_idx = row_ind[i] // min_size
    data_idx = row_ind[i] % min_size
    cluster_idx = col_ind[i]
    final_clusters[cluster_idx].append(datasets[dataset_idx][data_idx])

# Print the final clusters
for i, cluster in enumerate(final_clusters):
    print(f"Cluster {i + 1}:")
    for point in cluster:
        print(point)
