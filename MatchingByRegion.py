import time
import scanpy as sc
from scipy.sparse import lil_matrix
from scipy.spatial import KDTree
import numpy as np
import math
import sys
from lapjv import lapjv
import matplotlib.pyplot as plt


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


def weighted_distance(distance, scale, max_distance):
    """
    Applies an exponential weighting-up to the input distance, asymptotically approaching the max_distance.
    Scale adjusts the rate of growth.

    Higher scale more gradually increases the weighted distance
    Lower scale has less effect on lower distances, more effect on higher distances
    """

    # Apply exponential growth function
    return max_distance * (1 - math.exp(-distance / (scale * max_distance)))


def euclidean_distance(point1, point2):
    return np.sqrt(np.sum((point1 - point2) ** 2))


def match(adata1, adata2):
    if len(adata1) == 1 or len(adata2) == 1: return
    print(f"Started calculations for region {adata1[0].obs['CCFname']}...")

    # Extract spatial coordinates
    coords1 = extract_coords(adata1)
    coords2 = extract_coords(adata2)

    # # Plot this region in Coronal view
    # fig, neighours_fig = plt.subplots(figsize=(15, 10), dpi=800)
    # # Plot all cells from both datasets
    # neighours_fig.scatter(coords1[:, 2], coords1[:, 1], c='blue', label='adata1', alpha=0.5, s=0.2, linewidths=0.3)
    # neighours_fig.scatter(coords2[:, 2], coords2[:, 1], c='red', label='adata2', alpha=0.5, s=0.2, linewidths=0.3)
    # # for cell_Idx, cell_coord in enumerate(coords1):
    # #     neighours_fig.annotate(adata1.obs['slice'][cell_Idx], (cell_coord[2], cell_coord[1]), fontsize=1)
    # # for cell_Idx, cell_coord in enumerate(coords2):
    # #     neighours_fig.annotate(adata2.obs['slice'][cell_Idx], (cell_coord[2], cell_coord[1]), fontsize=1)
    # neighours_fig.set_xlabel('Z Coordinate')
    # neighours_fig.set_ylabel('Y Coordinate')
    # plt.savefig('graph.png', bbox_inches='tight')

    # Plot all cells from both datasets
    neighours_fig.scatter(coords1[:, 2], coords1[:, 1], c='blue', label='adata1', alpha=0.5, s=0.2, linewidths=0.2)
    neighours_fig.scatter(coords2[:, 2], coords2[:, 1], c='red', label='adata2', alpha=0.5, s=0.2, linewidths=0.2)
    label_cells = False
    if label_cells:
        for cell_Idx, cell_coord in enumerate(coords1):
            neighours_fig.annotate(cell_Idx, (cell_coord[0], cell_coord[1]), fontsize=4)

    # Construct KDTree for AnnData objects
    kdtree2 = KDTree(coords2)
    # Find neighboring cells
    num_valid_neighbours = num_neighbours
    if num_neighbours > len(coords2): num_valid_neighbours = len(coords2)
    distance_to_neighbours, neighbour_indices = kdtree2.query(coords1, k=num_valid_neighbours)
    # Find distance_threshold for every point in coords1
    distance_thresholds = distance_to_neighbours[:, -1]

    # Create matrix of gene expression x cells
    expression_matrix1 = adata1.X.toarray()
    expression_matrix2 = adata2.X.toarray()

    # Create matrix to hold all distances to each pair in adata2
    distances_matrix = np.full((len(coords1), len(coords2)), 2147483647, dtype=np.uint32)

    print("Calculating Distances...")
    time2 = time.time()
    for cell_idx in range(0, len(coords1)):
        neighbours_physical_matrix = coords2[neighbour_indices[cell_idx]]
        neighbours_expression_matrix = expression_matrix2[neighbour_indices[cell_idx]]

        physical_vector_norm = np.sum(coords1[cell_idx] ** 2)
        physical_norm2 = np.sum(neighbours_physical_matrix ** 2, axis=1)
        physical_dot_product = np.dot(neighbours_physical_matrix, coords1[cell_idx])
        physical_distances = np.sqrt(physical_vector_norm + physical_norm2 - 2 * physical_dot_product)

        expression_vector_norm = np.sum(expression_matrix1[cell_idx] ** 2)
        expression_norm2 = np.sum(neighbours_expression_matrix ** 2, axis=1)
        expression_dot_product = np.dot(neighbours_expression_matrix, expression_matrix1[cell_idx])
        expression_distances = np.sqrt(expression_vector_norm + expression_norm2 - 2 * expression_dot_product)

        for idx in range(0, len(neighbour_indices[cell_idx])):
            spatial_distance = physical_distances[idx]
            expression_distance = expression_distances[idx]
            if spatial_distance > 0 and expression_distance > 0:
                adjusted_spt_distance = weighted_distance(spatial_distance, 1, (distance_thresholds[cell_idx])) * 2
                adjusted_expression_distance = expression_distance * 100
                distances_matrix[cell_idx, neighbour_indices[cell_idx][idx]] = np.add(adjusted_spt_distance,
                                                                                      adjusted_expression_distance)

                # # Print results for debugging
                # print(f"\nCell {cell_idx} from adata1 to Cell {neighbour_indices[cell_idx][idx]} from adata2:\n"
                #     f"Calced Distance: {distances_matrix[cell_idx, neighbour_indices[cell_idx][idx]]}\n"
                #     # f"Exp Coords1: {','.join(map(str, expression_matrix1[cell_idx]))}\n"
                #     # f"Exp Coords2: {','.join(map(str, expression_matrix2[neighbour_indices[cell_idx, idx]]))}\n"
                #     f"Calced Exp Distance: {expression_distances[idx]}\n"
                #     f"Adjusted Exp Distance: {adjusted_expression_distance}\n"
                #     f"Spt Coords1: {','.join(map(str, coords1[cell_idx]))}\n"
                #     f"Spt Coords2: {','.join(map(str, coords2[neighbour_indices[cell_idx][idx]]))}\n"
                #     f"Calced Spt Distance: {physical_distances[idx]}\n"
                #     f"Adjusted Spt Distance: {adjusted_spt_distance}\n")

    print(f"Time to Calculate Euclidian Distances: {time.time() - time2}s")

    time3 = time.time()
    print("Finding best matches...")

    equalized_matrix = equalize_cardinality(distances_matrix)
    matches, col_index, _ = lapjv(equalized_matrix, verbose=0)

    global removed_cell_tracker
    removed_cells_before_region = removed_cell_tracker
    global all_matches
    global all_matches_tracker

    for idx in range(0, len(adata1)):
        if matches[idx] < len(adata2):
            cell_coords1 = coords1[idx, :]
            cell_coords2 = coords2[matches[idx], :]
            if valid_match(distances_matrix[idx, matches[idx]],
                           cell_coords1, cell_coords2, distance_thresholds[idx]):
                neighours_fig.plot([cell_coords1[2], cell_coords2[2]], [cell_coords1[1], cell_coords2[1]], 'r-', lw=0.03)
                all_matches[all_matches_tracker] = (adata1.obs_names[idx],
                                                    adata2.obs_names[matches[idx]])
                all_matches_tracker += 1
                # print("   Plotted cell.")
            else:
                removed_cell_tracker += 1
                # print("   Removed cell.")
        else:
            removed_cell_tracker += 1
            # print("   Removed cell.")

    print(f"Removed cells in this region: {removed_cell_tracker - removed_cells_before_region}")

    plt.savefig(plot_path, bbox_inches='tight')
    np.save('matches.npy', all_matches)

    print(f" Time to Calculate matches and plot: {time.time() - time3}s")
    print(f"Finished calculations for region {adata1[0].obs['CCFname']}\n\n")


def equalize_cardinality(matrix, fill_value=2147483647):
    """
    Pads a rectangular matrix to make it square by filling extra values with fill_value.

    Args:
    - matrix (list of lists or np.ndarray): The input matrix.
    - fill_value (int): The value to fill the extra cells. Default is 2147483647.

    Returns:
    - np.ndarray: A squared matrix.
    """
    if isinstance(matrix, list):
        matrix = np.array(matrix)

    rows, cols = matrix.shape
    size = max(rows, cols)

    # Create a new square matrix filled with the fill_value
    square_matrix = np.full((size, size), fill_value, dtype=matrix.dtype)

    # Copy the original matrix into the new square matrix
    square_matrix[:rows, :cols] = matrix

    return square_matrix


def valid_match(adjusted_distance, cell_coords1, cell_coords2, threshold):
    if adjusted_distance == 2147483647: return False
    distance = euclidean_distance(cell_coords1, cell_coords2)
    if distance <= threshold * 1.05 and distance <= absolute_distance_threshold: return True
    return False


start_time = time.time()
np.set_printoptions(edgeitems=100)

# Define constants
absolute_distance_threshold = 100  # In CCF units
num_neighbours = 300
plot_path = f"Plots/{time.time()}-AllMatches.png"

# Load AnnData objects
adata1_path = sys.argv[1]
adata2_path = sys.argv[2]

full_adata1 = sc.read_h5ad(adata1_path)
full_adata2 = sc.read_h5ad(adata2_path)

# Cut data into pieces for faster prototyping
cut_data = False
if cut_data:
    cut_data_factor = 40
    num_cells = full_adata1.shape[0]
    indices = np.random.permutation(num_cells)
    split = indices[:num_cells // cut_data_factor]
    full_adata1 = full_adata1[split].copy()

    num_cells = full_adata2.shape[0]
    indices = np.random.permutation(num_cells)
    split = indices[:num_cells // cut_data_factor]
    full_adata2 = full_adata2[split].copy()

    print('WARNING: Data split')

print(full_adata1)
print(full_adata2)

# # Plot the whole brain in 3d
# fig = plt.figure(figsize=(15, 10), dpi=900)
# neighours_fig = fig.add_subplot(111, projection='3d')
# # Plot all cells from both datasets
# coords1 = extract_coords((full_adata1))
# coords2 = extract_coords((full_adata2))
# neighours_fig.scatter(coords1[:, 0], coords1[:, 1], coords1[:, 2], c='blue', label='adata1', alpha=0.5, s=0.05,linewidths=0.05)
# neighours_fig.scatter(coords2[:, 0], coords2[:, 1], coords2[:, 2], c='red', label='adata2', alpha=0.5, s=0.05,linewidths=0.05)
# neighours_fig.set_xlabel('X Coordinate')
# neighours_fig.set_ylabel('Y Coordinate')
# neighours_fig.set_zlabel('Z Coordinate')
# neighours_fig.legend()
# plt.savefig(f'{time.time()}_graph.png', bbox_inches='tight')
# print("Big Plotting done")

# Create plot
fig, neighours_fig = plt.subplots(figsize=(15, 10), dpi=1200)

brain_regions1 = full_adata1.obs.groupby("CCFname")
brain_regions2 = full_adata2.obs.groupby("CCFname")

common_categories = sorted(set(brain_regions1.groups.keys()).intersection(set(brain_regions2.groups.keys())))

removed_cell_tracker = 0
all_matches = np.zeros(len(full_adata1), dtype=(np.int32, 2))
all_matches_tracker = 0

for category in common_categories:
    indices1 = brain_regions1.groups[category]
    indices2 = brain_regions2.groups[category]

    # adata1_subset = full_adata1[full_adata1.obs_names.isin(indices1)]
    # adata2_subset = full_adata2[full_adata2.obs_names.isin(indices2)]

    adata1_subset = full_adata1[indices1]
    adata2_subset = full_adata2[indices2]

    print(f"Calculating matches for region: {category}")
    print(adata1_subset)
    print(adata2_subset)
    match(adata1_subset, adata2_subset)

print(f"Total removed cells: {removed_cell_tracker}")
print(all_matches)

# Write matches to disk
np.save('matches.npy', all_matches)

neighours_fig.set_xlabel('X Coordinate')
neighours_fig.set_ylabel('Y Coordinate')
neighours_fig.legend()
print("Plotting Done.")

# ################## UMAP Matching #################
# # UMAPs are already computed in main.py
# # Get UMAP coords
# umap1 = adata1.obsm['X_umap']
# umap2 = adata2.obsm['X_umap']
#
# # Plot the UMAP embeddings
# fig, umap_fig = plt.subplots(figsize=(10, 10), dpi=500)
#
# # Plot UMAP for adata1
# umap_fig.scatter(umap1[:, 0], umap1[:, 1], c='blue', label='adata1', alpha=0.5, s=1)
# # Plot UMAP for adata2
# umap_fig.scatter(umap2[:, 0], umap2[:, 1], c='red', label='adata2', alpha=0.5, s=1)
#
# # Draw lines between matched cells
# for idx in range(len(adata1_match_idx)):
#     umap_coords1 = umap1[adata1_match_idx[idx], :]
#     umap_coords2 = umap2[adata2_match_idx[idx], :]
#     umap_fig.plot([umap_coords1[0], umap_coords2[0]], [umap_coords1[1], umap_coords2[1]], 'r-', lw=0.05)
#
# # Add labels and legend
# umap_fig.set_xlabel('UMAP1')
# umap_fig.set_ylabel('UMAP2')
# umap_fig.legend()
#
# plt.savefig('umap.png', bbox_inches='tight')

print(f'Script took: {time.time() - start_time}')