import time
import numpy
import scanpy as sc
import logging
from scipy.sparse import lil_matrix
from scipy.spatial import KDTree
import numpy as np
import math
import sys
import sslap

import PlotFromMatches


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

def weighted_distance(distance, scale):
    """
    Applies an exponential weighting-up to the input distance, asymptotically approaching the max_distance.
    Scale adjusts the rate of growth.

    Higher scale more gradually increases the weighted distance
    Lower scale has less effect on lower distances, more effect on higher distances
    """

    # Apply exponential growth function
    return distance_threshold * (1 - math.exp(-distance / (scale * distance_threshold)))

def match(adata1, adata2):

    # Extract spatial coordinates
    coords1 = extract_coords(adata1)
    coords2 = extract_coords(adata2)
    # coords1 = adata1.obsm['X_spatial']
    # coords2 = adata2.obsm['X_spatial']

    # Construct KDTree for both AnnData objects
    kdtree1 = KDTree(coords1)
    kdtree2 = KDTree(coords2)

    # Use query_ball_tree to find neighboring cells
    neighbour_indices = kdtree1.query_ball_tree(kdtree2, distance_threshold)

    # Find average number of neighbours
    total_items = sum(len(lst) for lst in neighbour_indices)
    number_of_lists = len(neighbour_indices)
    print(f"Average number of neighbours: {total_items / number_of_lists}")

    # Create matrix of gene expression x cells
    expression_matrix1 = adata1.X.toarray()
    expression_matrix2 = adata2.X.toarray()
    # expression_matrix1 = adata2.X
    # expression_matrix2 = adata2.X

    # Create matrix to hold all distances to each pair in adata2
    distances_matrix = lil_matrix((coords1.size, coords2.size))

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
                    weighted_distance(spatial_distance, 1.5) * 10, expression_distance / 1))

            # # Print results for debugging
            # print(f"\nCell {cell_idx} from adata1 to Cell {neighbour_indices[cell_idx][idx]} from adata2:\n"
            #       f"Calced Distance: {distances_matrix[cell_idx, neighbour_indices[cell_idx][idx]]}\n"
            #       # f"Exp Coords1: {','.join(map(str, expression_matrix1[cell_idx]))}\n"
            #       # f"Exp Coords2: {','.join(map(str, expression_matrix2[neighbour_indices[cell_idx, idx]]))}\n"
            #       f"Calced Exp Distance: {expression_distances[idx]}\n"
            #       f"Spt Coords1: {','.join(map(str, coords1[cell_idx]))}\n"
            #       f"Spt Coords2: {','.join(map(str, coords2[neighbour_indices[cell_idx][idx]]))}\n"
            #       f"Calced Spt Distance: {physical_distances[idx]}\n")

    print(f" Time to Calculate Euclidian Distances: {time.time() - time2}s")

    ################# Linear sum assignment method (no multiple mapping allowed) ################
    logger.info("Finding best matches...")
    coo_distances_matrix = distances_matrix.tocoo() * -1
    match_struct = sslap.auction_solve(coo_mat=coo_distances_matrix, problem='min', cardinality_check=False)
    matches = match_struct["sol"]

    # logger.info("Linear Sum Assignment solution:")
    # print(matches)

    return matches


start_time = time.time()

logger = logging.getLogger('logger')
grey = "\x1b[38;20m"
logging.basicConfig(level=logging.INFO, format="\x1b[38;20m" + "%(message)s")
np.set_printoptions(edgeitems=20)

# Define constants
distance_threshold = 3.5  # In CCF units

# Load AnnData objects
adata1_path = sys.argv[1]
adata2_path = sys.argv[2]

full_adata1 = sc.read_h5ad(adata1_path)
full_adata2 = sc.read_h5ad(adata2_path)

logger.info(full_adata1)
logger.info(full_adata2)

brain_regions1 = full_adata1.obs.groupby("CCFname")
brain_regions2 = full_adata2.obs.groupby("CCFname")

common_ccfnames = set(brain_regions1.groups.keys()).intersection(set(brain_regions2.groups.keys()))

matches = []

# Iterate over common CCFname values and apply do_work function to matched groups
for ccfname in common_ccfnames:
    matches += match(brain_regions1.get_group(ccfname), brain_regions2.get_group(ccfname))

unmatched_cells = 0
for idx in range(len(matches)):
    if not matches[idx] == -1:
        logger.info(f"{idx} {matches[idx]}")
    else:
        logger.info(f"REMOVED: {idx} {matches[idx]}")
        matches.pop(idx)
        unmatched_cells += 1
print(f"Unmatched Cells: {unmatched_cells}")

# # Write matches to disk
np.save('matches.npy', matches)

PlotFromMatches.plot(full_adata1, full_adata2, matches, distance_threshold)

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