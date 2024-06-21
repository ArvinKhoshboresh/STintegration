import time
import os
import scanpy as sc
import logging
from scipy.spatial import KDTree
from scipy.sparse.csgraph import min_weight_full_bipartite_matching
from scipy.sparse import csr_matrix, csc_matrix, lil_matrix
import matplotlib.pyplot as plt
import numpy as np
import math
import sys

start_time = time.time()

logger = logging.getLogger('logger')
logging.basicConfig(level=logging.INFO, format='')
np.set_printoptions(edgeitems=20)

# Load AnnData objects
adata1_path = sys.argv[1]
adata2_path = sys.argv[2]

adata1 = sc.read_h5ad(adata1_path)
adata2 = sc.read_h5ad(adata2_path)

logger.info(adata1)
logger.info(adata2)

# Check if spatial coordinates are present in both AnnData objects
if 'X_spatial' not in adata1.obsm:
    raise ValueError("Spatial coordinates not found in adata1.obsm['X_spatial']")
if 'X_spatial' not in adata2.obsm:
    raise ValueError("Spatial coordinates not found in adata2.obsm['X_spatial']")

###################### query ball point approach (in theory slower) ##################
# # Extract spatial coordinates and biuld KD tree
# coords1 = adata1.obsm['X_spatial']
# coords2 = adata2.obsm['X_spatial']
# kdtree = KDTree(coords2)
#
# # Define the distance threshold
# distance_threshold = 0.5
#
# # Iterate through all cells in the first AnnData object
# for idx, cell_id1 in enumerate(adata1.obs_names):
#     # Get the coordinates of the current cell
#     cell_coords1 = coords1[idx, :]
#
#     # Find all cells within the distance threshold in the second AnnData object
#     indices = kdtree.query_ball_point(cell_coords1, distance_threshold)
#
#     # Get the IDs of the cells within the threshold
#     cell_ids_within_threshold = adata2.obs_names[indices]
#
#     # print the result for the current cell
#     logger.info(f"Idx: {idx}, Cell ID from adata1: {cell_id1}")
#     logger.info(f"Cells within {distance_threshold} units in adata2: {cell_ids_within_threshold}")


######################### query ball tree approach (in theory faster) ##############################
# Extract spatial coordinates
coords1 = adata1.obsm['X_spatial']
coords2 = adata2.obsm['X_spatial']

# Construct KDTree for both AnnData objects
kdtree1 = KDTree(coords1)
kdtree2 = KDTree(coords2)

# Define constants
distance_threshold = 0.2  # In CCF units
available_memory = 60000000000  # In bytes, gets padded by 20%. 60GB
available_memory *= 0.8

# Use query_ball_tree to find neighboring cells
neighbour_indices = kdtree1.query_ball_tree(kdtree2, distance_threshold)

# Create plot
fig, neighours_fig = plt.subplots(figsize=(15, 10), dpi=500)

# Plot all cells from both datasets
neighours_fig.scatter(coords1[:, 0], coords1[:, 1], c='blue', label='adata1', alpha=0.5, s=1)
neighours_fig.scatter(coords2[:, 0], coords2[:, 1], c='red', label='adata2', alpha=0.5, s=1)
# for cell_Idx, cell_coord in enumerate(coords1):
#     neighours_fig.annotate(cell_Idx, (cell_coord[0], cell_coord[1]), fontsize=4)

# Create matrix of gene expression x cells
expression_matrix1 = adata1.X
expression_matrix2 = adata2.X

# Create matrix to hold all distances to each pair in adata2
distances_matrix = lil_matrix((coords1.shape[0], coords2.shape[0]))


def weighted_distance(distance, scale):
    """
    Applies an exponential weighting-up to the input distance, asymptotically approaching the max_distance.
    Scale adjusts the rate of growth.

    Higher scale more gradually increases the weighted distance
    Lower scale has less effect on lower distances, more effect on higher distances
    """

    # Apply exponential growth function
    return distance_threshold * (1 - math.exp(-distance / (scale * distance_threshold)))


# Equalize matrix shape
max_rows = min(expression_matrix1.shape[0], expression_matrix2.shape[0])
max_physical_columns = 3
max_expression_columns = min(expression_matrix1.shape[1], expression_matrix2.shape[1])

sliced_coords1 = coords1[:max_rows, :]
sliced_coords2 = coords2[:max_rows, :]
sliced_expression1 = expression_matrix1[:max_rows, :]
sliced_expression2 = expression_matrix2[:max_rows, :]

# Calculate number of chunks
element_size = sys.getsizeof(np.int32(0))  # Arbitrary int32 to get size
row_size = element_size * coords2.shape[0]
number_of_chunks = math.ceil(row_size * max_rows / available_memory)
rows_per_chunk = math.ceil(max_rows / number_of_chunks)

# Calculate chunks of dot products and prune for needed distances
for chunk in range(0, int(number_of_chunks)):
    chunked_coords1 = sliced_coords1[chunk * rows_per_chunk:chunk + 1 * rows_per_chunk, :]
    chunked_coords2 = sliced_coords2[chunk * rows_per_chunk:chunk + 1 * rows_per_chunk, :]

    chunked_expression1 = sliced_expression1[chunk * rows_per_chunk:chunk+1 * rows_per_chunk, :]
    chunked_expression2 = sliced_expression2[chunk * rows_per_chunk:chunk+1 * rows_per_chunk, :]

    # time1 = time.time()
    # physical_difference_matrix = np.subtract(chunked_coords1[:, np.newaxis, :], chunked_coords2[np.newaxis, :, :]) ** 2
    # physical_rank_order = np.sqrt(np.sum(physical_difference_matrix, axis=2))
    # print(f" time1: {time.time() - time1}")

    time2 = time.time()
    physical_norm1 = np.sum(chunked_coords1 ** 2, axis=1).reshape(-1, 1)
    physical_norm2 = np.sum(chunked_coords2 ** 2, axis=1).reshape(1, -1)
    physical_dot_product = np.dot(chunked_coords1, chunked_coords2.T)
    physical_distances = np.sqrt(physical_norm1 + physical_norm2 - 2 * physical_dot_product)

    expression_norm1 = np.sum(chunked_expression1 ** 2, axis=1).reshape(-1, 1)
    expression_norm2 = np.sum(chunked_expression2 ** 2, axis=1).reshape(1, -1)
    expression_dot_product = np.dot(chunked_expression1, chunked_expression2.T)
    expression_distances = np.sqrt(expression_norm1 + expression_norm2 - 2 * expression_dot_product)
    print(f" Time to Calculate Euclidian Distances: {time.time() - time2}s")

    expression_difference_matrix = np.subtract(chunked_expression1, chunked_expression2) ** 2
    expression_rank_order = expression_difference_matrix.dot(expression_difference_matrix.T)

    neighbours_removed = 0
    for cell_idx in range(0, rows_per_chunk):
            cell_neighbours = neighbour_indices[(max_rows * chunk) + cell_idx]
            for neighbour_idx in cell_neighbours:
                try:
                    weighted_physical_distance = weighted_distance(physical_distances[cell_idx, neighbour_idx], 2)
                    expression_distance = expression_distances[cell_idx, neighbour_idx]
                    distances_matrix[(max_rows * chunk) + cell_idx, neighbour_idx] = np.add(weighted_physical_distance * 100000,
                                                                                            expression_distance / 100000)

                    # print(f"\nCell {(max_rows * chunk) + cell_idx} from adata1 to Cell {neighbour_idx} from adata2:\n"
                    #       f"Calced Distance: {distances_matrix[(max_rows * chunk) + cell_idx, neighbour_idx]}\n"
                    #       # f"Exp Coords1: {','.join(map(str, chunked_expression1[cell_idx]))}\n"
                    #       # f"Exp Coords2: {','.join(map(str, chunked_expression2[neighbour_idx]))}\n"
                    #       f"Calced Exp Distance: {expression_distances[cell_idx, neighbour_idx]}\n"
                    #       f"Spt Coords1: {','.join(map(str, chunked_coords1[cell_idx]))}\n"
                    #       f"Spt Coords2: {','.join(map(str, chunked_coords2[neighbour_idx]))}\n"
                    #       f"Calced Spt Distance: {physical_distances[cell_idx, neighbour_idx]}\n")

                except:
                    neighbours_removed += 1
    logger.info(f"Neighbours Removed when Equalizing Shape: {neighbours_removed}")

# ####################### FOR LOOP METHOD #########################
#
#
# def dot_product_euclidean_distance(coords1, coords2):
#     """
#     Takes two cell coords in many dimensional space, returns Euclidean distance approximated by dot product of
#     difference vector
#     """
#
#     difference_vectors = coords1 - coords2
#     return np.sqrt(np.matmul(difference_vectors, difference_vectors))
#
#
# for_loop_distances_matrix = lil_matrix((len(coords1), len(coords2)), dtype=np.int32)
#
# # Iterate through all cells in adata1 and plot neighbors from adata2
# for idx, cell_id1 in enumerate(adata1.obs_names):
#
#     # Get coordinates of the current cell
#     cell_coords1 = coords1[idx, :]
#
#     # Get indices of neighboring cells in adata2
#     cell_neighbours = neighbour_indices[idx]
#
#     # Get vector representation of adata1 cell's gene expression
#     expression_vector1 = expression_matrix1[idx]
#
#     # print result for the current cell
#     #logger.info(f"Idx: {idx}, Cell ID from adata1: {cell_id1}")
#     #logger.info(f"Cells within {distance_threshold} units in adata2: \n")
#
#     # Plot lines connecting the current cell to its neighbors and calculate distances
#     for neighbour_idx in cell_neighbours:
#         cell_coords2 = coords2[neighbour_idx, :]
#
#         # Calculate distance and weighted distance
#         weighted_physical_distance = weighted_distance(dot_product_euclidean_distance(cell_coords1, cell_coords2), 2)
#         weighted_physical_distance = weighted_physical_distance * 100000  # Creates values broadly in the 3 figures
#
#         # Add distance to physical distance matrix
#         for_loop_distances_matrix[idx, neighbour_idx] = weighted_physical_distance
#
#         # Plot match line
#         # neighours_fig.plot([cell_coords1[0], cell_coords2[0]], [cell_coords1[1], cell_coords2[1]], 'k-', lw=0.25)
#
#         # print distance for verification
#         #logger.info(f"Distance from {cell_id1} to {adata2.obs_names[neighbour_idx]}: {weighted_physical_distance:.6f}")
#
#         # Get match expression vector
#         expression_vector2 = expression_matrix2[neighbour_idx]
#
#         expression_distance = dot_product_euclidean_distance(expression_vector1, expression_vector2)
#         for_loop_distances_matrix[idx, neighbour_idx] = ((np.add(distances_matrix[idx, neighbour_idx],
#                                                                  expression_distance)) / 100000)
#
# for idx in range(0, coords1.shape[0]):
#     for idx2 in range(0, coords2.shape[0]):
#         first = distances_matrix[idx, idx2]
#         second = for_loop_distances_matrix[idx, idx2]
#         if first or second:
#             print(f"Matrix Method: {first}   For-Loop Method: {second}")
#
###################### END FOR LOOP METHOD (GOOD RIDDANCE) ######################


distances_matrix *= -1
distances_matrix = distances_matrix.tocsr()

# ############### MULTIPLE MAPPING ALLOWED METHOD ################
# # Draw red line for lowest expression Euclidean distance match
# for idx, cell_id1 in enumerate(adata1.obs_names):
#     # Get index of the closest neighbor
#     closest_neighbor_idx = np.argmin(physical_distances[idx, :])
#     if physical_distances[idx, closest_neighbor_idx] != 2147483647:
#         cell_coords1 = coords1[idx, :]
#         cell_coords2 = coords2[closest_neighbor_idx, :]
#
#         # Plot red line for the closest neighbor
#         neighours_fig.plot([cell_coords1[0], cell_coords2[0]], [cell_coords1[1], cell_coords2[1]], 'r-', lw=0.5)
# # TODO: Add save to disk functionality

################# Linear sum assingment method (no multiple mapping allowed) ################
logger.info("Finding best matches...")
adata1_match_idx, adata2_match_idx = min_weight_full_bipartite_matching(distances_matrix)
logger.info("Linear Sum Assignment solution:")
for idx in range(len(adata1_match_idx)):
    if dot_product_euclidean_distance(coords1[adata1_match_idx[idx]],
                                      coords2[adata2_match_idx[idx]]) <= distance_threshold:
        cell_coords1 = coords1[adata1_match_idx[idx], :]
        cell_coords2 = coords2[adata2_match_idx[idx], :]
        neighours_fig.plot([cell_coords1[0], cell_coords2[0]], [cell_coords1[1], cell_coords2[1]], 'r-', lw=0.05)
        logger.info(f"{adata1_match_idx[idx]} {adata2_match_idx[idx]}")
    else:
        adata1_match_idx.pop(idx)
        adata2_match_idx.pop(idx)

# Write matches to disk
matches_array_length = max(len(adata1_match_idx), len(adata2_match_idx))
with open('matches.txt', 'w') as file:
    for idx in range(matches_array_length):
        elem1 = adata1_match_idx[idx] if idx < len(adata1_match_idx) else None
        elem2 = adata2_match_idx[idx] if idx < len(adata2_match_idx) else None
        file.write(f"{elem1} {elem2}\n")

# Add labels and legend
neighours_fig.set_xlabel('X Coordinate')
neighours_fig.set_ylabel('Y Coordinate')
neighours_fig.legend()

plt.savefig('matches.png', bbox_inches='tight')

################## UMAP Matching #################
# UMAPs are already computed in main.py
# Get UMAP coords
umap1 = adata1.obsm['X_umap']
umap2 = adata2.obsm['X_umap']

# Plot the UMAP embeddings
fig, umap_fig = plt.subplots(figsize=(10, 10), dpi=500)

# Plot UMAP for adata1
umap_fig.scatter(umap1[:, 0], umap1[:, 1], c='blue', label='adata1', alpha=0.5, s=1)
# Plot UMAP for adata2
umap_fig.scatter(umap2[:, 0], umap2[:, 1], c='red', label='adata2', alpha=0.5, s=1)

# Draw lines between matched cells
for idx in range(len(adata1_match_idx)):
    umap_coords1 = umap1[adata1_match_idx[idx], :]
    umap_coords2 = umap2[adata2_match_idx[idx], :]
    umap_fig.plot([umap_coords1[0], umap_coords2[0]], [umap_coords1[1], umap_coords2[1]], 'r-', lw=0.05)

# Add labels and legend
umap_fig.set_xlabel('UMAP1')
umap_fig.set_ylabel('UMAP2')
umap_fig.legend()

plt.savefig('umap.png', bbox_inches='tight')

logger.info(f'Script took: {time.time() - start_time}')
