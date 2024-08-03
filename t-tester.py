import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
import anndata


# Placeholder function for extracting coordinates from AnnData
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


# Sample data
# Replace this with your actual list of lists, with each sublist containing 8 elements
list_of_lists = [
    [1, 2, 3, 4, 5, 6, 7, 8],
    [9, 10, 11, 12, 13, 14, 15, 16]
]

# List of your actual AnnData file paths
adata_files = [
    'file1.h5ad', 'file2.h5ad', 'file3.h5ad', 'file4.h5ad',
    'file5.h5ad', 'file6.h5ad', 'file7.h5ad', 'file8.h5ad'
]

# Load all AnnData objects
adatas = [anndata.read_h5ad(file) for file in adata_files]

# Extract coordinates from all AnnData objects
all_coords = [extract_coords(adata) for adata in adatas]

t_stats = []
average_coords = []

for sublist in list_of_lists:
    group1_ids = sublist[:4]
    group2_ids = sublist[4:]

    group1_expr = [adatas[i][adatas[i].obs.index == str(id)].X.flatten() for i, id in enumerate(group1_ids)]
    group2_expr = [adatas[i][adatas[i].obs.index == str(id)].X.flatten() for i, id in enumerate(group2_ids)]

    # Perform t-test
    t_stat, p_value = ttest_ind(group1_expr, group2_expr)
    t_stats.append(t_stat)

    # Get average coordinates
    coords = [all_coords[i][id] for i, id in enumerate(sublist)]
    avg_coords = np.mean(coords, axis=0)
    average_coords.append(avg_coords)

# Normalize t-stats for color mapping
norm = plt.Normalize(min(t_stats), max(t_stats))
cmap = plt.get_cmap('coolwarm')

# Plotting
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for avg_coord, t_stat in zip(average_coords, t_stats):
    color = cmap(norm(t_stat))
    ax.scatter(avg_coord[0], avg_coord[1], avg_coord[2], color=color)

# Add a color bar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
plt.colorbar(sm, ax=ax, label='t-statistic')

plt.show()
