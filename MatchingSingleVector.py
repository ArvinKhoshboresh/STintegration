import time
import numpy
import scanpy as sc
import logging
from scipy.sparse import lil_matrix
from scipy.spatial import KDTree
import matplotlib.pyplot as plt
import numpy as np
import math
import sys
import sslap

start_time = time.time()

logger = logging.getLogger('logger')
grey = "\x1b[38;20m"
logging.basicConfig(level=logging.INFO, format="\x1b[38;20m" + "%(message)s")
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

# Use query_ball_tree to find neighboring cells
neighbour_indices = kdtree1.query_ball_tree(kdtree2, distance_threshold)

# Create plot
fig, neighours_fig = plt.subplots(figsize=(15, 10), dpi=800)

# Plot all cells from both datasets
neighours_fig.scatter(coords1[:, 0], coords1[:, 1], c='blue', label='adata1', alpha=0.5, s=0.2, linewidths=0.3)
neighours_fig.scatter(coords2[:, 0], coords2[:, 1], c='red', label='adata2', alpha=0.5, s=0.2, linewidths=0.3)
# for cell_Idx, cell_coord in enumerate(coords1):
#     neighours_fig.annotate(cell_Idx, (cell_coord[0], cell_coord[1]), fontsize=4)

# Create matrix of gene expression x cells
expression_matrix1 = adata1.X
expression_matrix2 = adata2.X

# Create matrix to hold all distances to each pair in adata2
# distances_matrix = np.full((coords1.shape[0], coords2.shape[0]), fill_value=32767, dtype=np.int16)
distances_matrix = lil_matrix((coords1.size, coords2.size))


def weighted_distance(distance, scale):
    """
    Applies an exponential weighting-up to the input distance, asymptotically approaching the max_distance.
    Scale adjusts the rate of growth.

    Higher scale more gradually increases the weighted distance
    Lower scale has less effect on lower distances, more effect on higher distances
    """

    # Apply exponential growth function
    return distance_threshold * (1 - math.exp(-distance / (scale * distance_threshold)))

print("Calculating Distances...")
time2 = time.time()
for cell_idx in range(0, len(coords1)):
    neighbours_physical_matrix = numpy.zeros((len(neighbour_indices[cell_idx]), len(coords2[0])))
    neighbours_expression_matrix = numpy.zeros((len(neighbour_indices[cell_idx]), len(expression_matrix2[0])))
    for idx, neighbour in enumerate(neighbour_indices[cell_idx]):
        neighbours_physical_matrix[idx] = coords2[neighbour]
        neighbours_expression_matrix[idx] = expression_matrix2[neighbour]

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
            distances_matrix[cell_idx, neighbour_indices[cell_idx][idx]] = round(np.add(
                weighted_distance(spatial_distance, 1.5) * 1000, expression_distance / 30000))

        # print(f"\nCell {cell_idx} from adata1 to Cell {neighbour_indices[cell_idx][idx]} from adata2:\n"
        #       f"Calced Distance: {distances_matrix[cell_idx, neighbour_indices[cell_idx][idx]]}\n"
        #       # f"Exp Coords1: {','.join(map(str, expression_matrix1[cell_idx]))}\n"
        #       # f"Exp Coords2: {','.join(map(str, expression_matrix2[neighbour_indices[cell_idx, idx]]))}\n"
        #       f"Calced Exp Distance: {expression_distances[idx]}\n"
        #       f"Spt Coords1: {','.join(map(str, coords1[cell_idx]))}\n"
        #       f"Spt Coords2: {','.join(map(str, coords2[neighbour_indices[cell_idx][idx]]))}\n"
        #       f"Calced Spt Distance: {physical_distances[idx]}\n")

print(f" Time to Calculate Euclidian Distances: {time.time() - time2}s")

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
coo_distances_matrix = distances_matrix.tocoo() * -1
match_struct = sslap.auction_solve(coo_mat=coo_distances_matrix, problem='min', cardinality_check=False)
matches = match_struct["sol"]
logger.info("Linear Sum Assignment solution:")

unmatched_cells = 0

for idx in range(len(matches)):
    if not matches[idx] == -1:
        logger.info(f"{idx} {matches[idx]}")
        cell_coords1 = coords1[idx, :]
        cell_coords2 = coords2[matches[idx], :]
        neighours_fig.plot([cell_coords1[0], cell_coords2[0]], [cell_coords1[1], cell_coords2[1]], 'r-', lw=0.03)
    else:
        unmatched_cells += 1


# # Write matches to disk
print(f"Unmatched Cells: {unmatched_cells}")
np.save('matches.npy', matches)
# matches_array_length = max(len(adata1_match_idx), len(adata2_match_idx))
# with open('matches.txt', 'w') as file:
#     for idx in range(matches_array_length):
#         elem1 = adata1_match_idx[idx] if idx < len(adata1_match_idx) else None
#         elem2 = adata2_match_idx[idx] if idx < len(adata2_match_idx) else None
#         file.write(f"{elem1} {elem2}\n")

# Add labels and legend
neighours_fig.set_xlabel('X Coordinate')
neighours_fig.set_ylabel('Y Coordinate')
neighours_fig.legend()

plt.savefig('matches.png', bbox_inches='tight')

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

logger.info(f'Script took: {time.time() - start_time}')
