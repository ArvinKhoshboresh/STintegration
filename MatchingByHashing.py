import numpy as np
import sys
import scanpy as sc
import numpy as np
import falconn
from typing import List, Tuple


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
    cut_data_factor = 10000
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


def lsh_vector_groups(vector_lists: List[List[np.ndarray]], num_groups: int = 5) -> List[Tuple[np.ndarray]]:
    # Combine all vectors and keep track of their origins
    all_vectors = []
    vector_origins = []
    for list_idx, vector_list in enumerate(vector_lists):
        all_vectors.extend(vector_list)
        vector_origins.extend([list_idx] * len(vector_list))

    all_vectors = np.array(all_vectors)

    # Set up LSH parameters
    params_cp = falconn.LSHConstructionParameters()
    params_cp.dimension = all_vectors.shape[1]
    params_cp.lsh_family = falconn.LSHFamily.CrossPolytope
    params_cp.distance_function = falconn.DistanceFunction.EuclideanSquared
    params_cp.l = 10  # Number of hash tables
    params_cp.num_rotations = 2
    params_cp.seed = 5721840
    params_cp.num_setup_threads = 0
    params_cp.k = 1
    params_cp.storage_hash_table = falconn.StorageHashTable.BitPackedFlatHashTable
    params_cp.last_cp_dimension = 1

    # Set up LSH index
    table = falconn.LSHIndex(params_cp)
    table.setup(all_vectors)

    query_object = table.construct_query_object()

    # Find groups
    groups = []
    for i in range(len(all_vectors)):
        nearest = query_object.get_unique_candidates(all_vectors[i])

        # Check if we have vectors from all lists
        if len(set(vector_origins[j] for j in nearest)) == len(vector_lists):
            # Form a group
            group = [all_vectors[j] for j in nearest if vector_origins[j] == list_idx][:1]
            for list_idx in range(len(vector_lists)):
                group.extend([all_vectors[j] for j in nearest if vector_origins[j] == list_idx][:1])

            if len(group) == len(vector_lists):
                groups.append(tuple(group))

        if len(groups) >= num_groups:
            break

    # Sort groups by total pairwise distance
    def total_distance(group):
        return sum(np.sum((v1 - v2) ** 2) for v1 in group for v2 in group)

    groups.sort(key=total_distance)

    return groups[:num_groups]


# Create some sample data
np.random.seed(42)
vector_lists = [
    [np.random.rand(10) for _ in range(100)],
    [np.random.rand(10) for _ in range(100)],
    [np.random.rand(10) for _ in range(100)]
]

result = lsh_vector_groups(vector_lists, num_groups=3)

print(result)

for i, group in enumerate(result):
    print(f"Group {i + 1}:")
    for vector in group:
        print(vector)
