import sys
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from line_profiler_pycharm import profile


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


def permute_rows(matrix):
    # Convert input matrix to numpy array if it isn't already
    matrix = np.array(matrix)

    # Get the number of rows and columns
    n, m = matrix.shape

    # Create a new matrix to store the randomized rows
    randomized_matrix = np.empty_like(matrix)

    # Randomly shuffle rows for each column
    for col in range(m):
        randomized_matrix[:, col] = np.random.permutation(matrix[:, col])

    return randomized_matrix


@profile
def main():
    np.set_printoptions(edgeitems=25)

    # Load AnnData objects
    adata_paths = [sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8]]
    adatas = [sc.read(file) for file in adata_paths]
    print(adatas)

    all_coords = [extract_coords(adata) for adata in adatas]

    match_chains = sys.argv[9]

    chains = np.load(match_chains)
    print(f"Matches length: {len(chains)}")

    cut_chains = True
    if cut_chains:
        np.random.seed(42)
        cut_data_factor = 20000
        num_chains = len(chains)
        indices = np.random.permutation(num_chains)
        split = indices[:num_chains // cut_data_factor]
        chains = chains[split].copy()
        print(f"WARNING: Chains cut to: {len(chains)}")

    permute_chains = False
    if permute_chains:
        chains = permute_rows(chains)

    print(chains)

    t_stats = []
    average_coords = []

    for chain in chains:
        control_ids = chain[:4]
        enucleated_ids = chain[4:]

        control_expr = [adatas[i][str(idx)].X.toarray() for i, idx in enumerate(control_ids)]
        enucleated_expr = [adatas[i+4][str(idx)].X.toarray() for i, idx in enumerate(enucleated_ids)]

        # Perform t-test
        t_stat, p_value = ttest_ind(control_expr, enucleated_expr)
        print(t_stat)
        t_stats.append(t_stat)

        # Get average coordinates
        coords = [all_coords[i][np.where(adatas[i].obs_names == str(idx))] for i, idx in enumerate(chain)]
        avg_coords = np.mean(coords, axis=0)
        print(f"avg_coords:{avg_coords}")
        average_coords.append(avg_coords)

    # Normalize t-stats for color mapping
    norm = plt.Normalize(min(t_stats), max(t_stats))
    cmap = plt.get_cmap('coolwarm')

    # Plotting
    fig = plt.figure(dpi=1500)
    ax = fig.add_subplot(111, projection='3d')

    for avg_coord, t_stat in zip(average_coords, t_stats):
        color = cmap(norm(t_stat))
        ax.scatter(avg_coord[0], avg_coord[1], avg_coord[2], color=color)

    # Add a color bar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, ax=ax, label='t-statistic')

    plt.savefig("warmth plot", bbox_inches='tight')


if __name__ == "__main__":
    main()
