import time
import os
import scanpy as sc
import logging
from scipy.spatial import KDTree
from scipy.optimize import linear_sum_assignment
import matplotlib.pyplot as plt
import numpy as np
import math
import sys

start_time = time.time()

logger = logging.getLogger('logger')
logging.basicConfig(level=logging.INFO, format='')

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

# Define distance threshold
distance_threshold = 0.3

# Use query_ball_tree to find neighboring cells
neighbour_indices = kdtree1.query_ball_tree(kdtree2, distance_threshold)

# Create matrix to hold all physical distances to each pair in adata2
physical_distances = np.full((coords1.shape[0], coords2.shape[0]), fill_value=0, dtype=np.uint16)

# Create matrix to hold all adata1 cells and expression Euclidean distance to each cell in adata2
expression_distances = np.full((coords1.shape[0], coords2.shape[0]), fill_value=0, dtype=np.uint16)

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


def dot_product_euclidean_distance(coords1, coords2):
    """
    Takes two cell coords in many dimensional space, returns Euclidean distance approximated by dot product of
    difference vector
    """

    difference_vectors = coords1 - coords2
    return np.sqrt(np.dot(difference_vectors, difference_vectors))


def weighted_distance(distance, scale):
    """
    Applies an exponential weighting-up to the input distance, asymptotically approaching the max_distance.
    Scale adjusts the rate of growth.

    Higher scale more gradually increases the weighted distance
    Lower scale has less effect on lower distances, more effect on higher distances
    """

    # Apply exponential growth function
    return distance_threshold * (1 - math.exp(-distance / (scale * distance_threshold)))


# Iterate through all cells in adata1 and plot neighbors from adata2
for idx, cell_id1 in enumerate(adata1.obs_names):

    # Get coordinates of the current cell
    cell_coords1 = coords1[idx, :]

    # Get indices of neighboring cells in adata2
    cell_neighbours = neighbour_indices[idx]

    # Get vector representation of adata1 cell's gene expression
    expression_vector1 = expression_matrix1[idx]

    # print result for the current cell
    logger.info(f"Idx: {idx}, Cell ID from adata1: {cell_id1}")
    logger.info(f"Cells within {distance_threshold} units in adata2: \n")

    # Plot lines connecting the current cell to its neighbors and calculate distances
    for neighbour_idx in cell_neighbours:
        cell_coords2 = coords2[neighbour_idx, :]

        # Calculate distance and weighted distance
        weighted_physical_distance = weighted_distance(dot_product_euclidean_distance(cell_coords1, cell_coords2), 2)
        weighted_physical_distance = weighted_physical_distance * 1000  # Creates values broadly in the 3 figures

        # Add distance to physical distance matrix
        physical_distances[idx, neighbour_idx] = weighted_physical_distance

        # Plot match line
        # neighours_fig.plot([cell_coords1[0], cell_coords2[0]], [cell_coords1[1], cell_coords2[1]], 'k-', lw=0.25)

        # print distance for verification
        logger.info(f"Distance from {cell_id1} to {adata2.obs_names[neighbour_idx]}: {weighted_physical_distance:.6f}")

        expression_vector2 = expression_matrix2[neighbour_idx]

        expression_distance = dot_product_euclidean_distance(expression_vector1, expression_vector2)
        expression_distances[idx, neighbour_idx] = expression_distance / 1000  # Creates values broadly in the 3 figures

np.set_printoptions(threshold=1000)

logger.info("Physical Matrix:")
logger.info(physical_distances)

logger.info("Expression Matrix:")
logger.info(expression_distances)

distances_matrix = np.add(physical_distances, expression_distances)
distances_matrix[distances_matrix == 0] = 65535

logger.info("Distances Matrix:")
logger.info(distances_matrix)


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

#################### Linear sum assingment method (no multiple mapping allowed) ################
adata1_match_idx, adata2_match_idx = linear_sum_assignment(distances_matrix)

logger.info("Linear Sum Assignment solution:")

for idx in range(len(adata1_match_idx)):
    if physical_distances[adata1_match_idx[idx], adata2_match_idx[idx]] <= distance_threshold:
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
