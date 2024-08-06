import sys
import time
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from statsmodels.multivariate.manova import MANOVA
from sklearn.decomposition import PCA
import concurrent.futures

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


def array_transformer(array):

    min_val = np.min(array)
    max_val = np.max(array)
    array = (array - min_val) / (max_val - min_val)

    array = np.log(array)

    return array


def manova_calculator(idx, chain):

    control_ids = chain[:4]
    enucleated_ids = chain[4:]

    controls = [adatas[idx2][obs_names_dict[idx2][str(control_id)]]
                for idx2, control_id in enumerate(control_ids)]
    enucleateds = [adatas[idx2 + 4][obs_names_dict[idx2 + 4][str(enucleated_id)]]
                   for idx2, enucleated_id in enumerate(enucleated_ids)]

    control_expr = [row.X.toarray()[0]
                    for row in controls]
    enucleated_expr = [row.X.toarray()[0]
                       for row in enucleateds]

    data = np.vstack([control_expr, enucleated_expr])

    n_components = min(data.shape[0], data.shape[1]) - 1
    pca = PCA(n_components=n_components)
    data = pca.fit_transform(data)

    labels = ['Control'] * len(control_expr) + ['Enucleated'] * len(enucleated_expr)
    gene_numbers = [f'Gene{i + 1}' for i in range(n_components)]
    df = pd.DataFrame(data, columns=gene_numbers)
    df['Group'] = labels

    formula = ' + '.join(gene_numbers) + ' ~ Group'
    manova = MANOVA.from_formula(formula, data=df)
    try:
        result = manova.mv_test()
    except:
        print(f"SKIPPED {idx}. Probably a Singular Matrix.")

    hotelling_lawley_trace = result.results['Group']['stat'].loc['Hotelling-Lawley trace', 'Value']

    # Get average coordinates
    control_coords = [extract_coords(row) for row in controls]
    enucleated_coords = [extract_coords(row) for row in enucleateds]
    coords = np.vstack(control_coords + enucleated_coords)
    avg_coords = np.mean(coords, axis=0)

    return idx, avg_coords, hotelling_lawley_trace


def main():
    time1 = time.time()
    np.set_printoptions(edgeitems=35)

    # Load AnnData objects
    adata_paths = [sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],
                   sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8]]
    global adatas
    adatas = [sc.read(file) for file in adata_paths]
    print(adatas)

    # all_coords = [extract_coords(adata) for adata in adatas]

    match_chains = sys.argv[9]

    chains = np.load(match_chains)
    global len_chains
    len_chains = len(chains)
    print(f"Matches length: {len_chains}")

    cut_chains = True
    if cut_chains:
        np.random.seed(42)
        cut_data_factor = 1000
        num_chains = len(chains)
        indices = np.random.permutation(num_chains)
        split = indices[:num_chains // cut_data_factor]
        chains = chains[split].copy()
        len_chains = len(chains)
        print(f"WARNING: Chains cut to: {len_chains}")

    permute_chains = False
    if permute_chains:
        chains = permute_rows(chains)

    print(chains)

    global obs_names_dict

    manova_stats = []
    average_coords = []

    obs_names_dict = [{name: idx for idx, name in enumerate(adata.obs_names)} for adata in adatas]

    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = [executor.submit(manova_calculator, idx, chain) for idx, chain in enumerate(chains)]

        for future in concurrent.futures.as_completed(futures):
            idx, avg_coords, manova_stat = future.result()
            print(f"{idx + 1}/{len_chains}")
            print(f"Avg. Coords: {avg_coords}")
            print(f"MANOVA Stat: {manova_stat}")
            average_coords.append(avg_coords)
            manova_stats.append(manova_stat)

    # norm = plt.Normalize(vmin=min(manova_stats), vmax=max(manova_stats))
    cmap = plt.get_cmap('coolwarm')
    print(np.array(manova_stats))
    manova_stats = np.array(manova_stats)
    max_non_inf = np.max(manova_stats[np.isfinite(manova_stats)])
    manova_stats[np.isinf(manova_stats)] = max_non_inf
    manova_stats = array_transformer(manova_stats)
    print(manova_stats)

    # TODO: COLORS STILL NOT WORKING

    # Plotting
    fig = plt.figure(dpi=1500)
    ax = fig.add_subplot(111, projection='3d')

    for avg_coord, manova_stat in zip(average_coords, manova_stats):
        color = cmap(manova_stat)
        ax.scatter(avg_coord[0], avg_coord[1], avg_coord[2], color=color, s=0.1, lw=0)
        # ax.text(avg_coord[0], avg_coord[1], avg_coord[2], manova_stat, fontsize=3)

    plt.savefig("warmth plot", bbox_inches='tight')

    print("Plotted.")
    print(f"Script took: {time.time() - time1}")


if __name__ == "__main__":
    main()

len_chains = 0
adatas = []
obs_names_dict = {}
