import numpy as np
import sys
import logging
import matplotlib.pyplot as plt
import scanpy as sc
import MatchingByRegion


def euclidean_distance(point1, point2):
    return np.sqrt(np.sum((point1 - point2) ** 2))


def plot(adata1, adata2, matches, threshold):

    print("Started plotting...")

    # Extract spatial coordinates
    coords1 = MatchingByRegion.extract_coords(adata1)
    coords2 = MatchingByRegion.extract_coords(adata2)
    # coords1 = adata1.obsm['X_spatial']
    # coords2 = adata2.obsm['X_spatial']

    # Create plot
    fig, neighours_fig = plt.subplots(figsize=(15, 10), dpi=800)

    # Plot all cells from both datasets
    neighours_fig.scatter(coords1[:, 0], coords1[:, 1], c='blue', label='adata1', alpha=0.5, s=0.2, linewidths=0.3)
    neighours_fig.scatter(coords2[:, 0], coords2[:, 1], c='red', label='adata2', alpha=0.5, s=0.2, linewidths=0.3)

    label_cells = False
    if label_cells:
        for cell_Idx, cell_coord in enumerate(coords1):
            neighours_fig.annotate(cell_Idx, (cell_coord[0], cell_coord[1]), fontsize=4)

    removed_cells = 0

    # Loop through matches and plot
    for idx in range(len(matches)):
        if matches[idx] != -1:
            logger.info(f"{idx} {matches[idx]}")
            cell_coords1 = coords1[idx, :]
            cell_coords2 = coords2[matches[idx], :]
            distance = euclidean_distance(cell_coords1, cell_coords2)
            if distance < threshold:
                neighours_fig.plot([cell_coords1[0], cell_coords2[0]], [cell_coords1[1], cell_coords2[1]], 'r-',
                                   lw=0.03)
            else:
                removed_cells += 1

    print(f"Removed Cells: {removed_cells}")

    neighours_fig.set_xlabel('X Coordinate')
    neighours_fig.set_ylabel('Y Coordinate')
    neighours_fig.legend()

    plt.savefig('matches.png', bbox_inches='tight')

    print("Plotting Done.")


# Setup logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()

# Load AnnData objects
adata1_path = sys.argv[1]
adata2_path = sys.argv[2]

full_adata1 = sc.read_h5ad(adata1_path)
full_adata2 = sc.read_h5ad(adata2_path)

logger.info(full_adata1)
logger.info(full_adata2)

distance_threshold = 3.5

# Load arrays from disk
loaded_matches = np.load('matches.npy')

plot(full_adata1, full_adata2, loaded_matches, distance_threshold)

print("Script Completed")
