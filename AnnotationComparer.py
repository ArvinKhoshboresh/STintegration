import time

import numpy as np
import scanpy as sc
import sys
import random
from collections import Counter


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


# Load AnnData objects
adata_paths = [sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8]]
adata_struct = [sc.read(file) for file in adata_paths]

match_chains = sys.argv[9]

np.set_printoptions(edgeitems=25)

print(adata_struct)

chains = np.load(match_chains)
print(f"Matches length: {len(chains)}")

cut_chains = True
if cut_chains:
    np.random.seed(42)
    cut_data_factor = 100
    num_chains = len(chains)
    indices = np.random.permutation(num_chains)
    split = indices[:num_chains // cut_data_factor]
    chains = chains[split].copy()
    print(f"WARNING: Chains cut to: {len(chains)}")

permute_chains = False
if permute_chains:
    chains = permute_rows(chains)

chains_length = len(chains)
print(chains)

annotation_names = ['clustid', 'clustname', 'subclass']

time1 = time.time()
# annotations_list = [[list(adata_struct[idx][str(obs_name)].obs[annotation_names].values[0])
#                     for idx, obs_name in enumerate(chain)]
#                     for chain in chains]
annotations_list = np.array([[adata_struct[idx][str(obs_name)].obs[annotation_names].values[0]
                              for idx, obs_name in enumerate(chain)]
                             for chain in chains])

# print(annotations_list)

total_averages = [0, 0, 0, 0, 0]
for idx, annotations in enumerate(annotations_list):
    annotations_length = len(annotations)
    for idx2 in range(len(total_averages)):
        if idx2 < 3:
            largest_matching_count = max(Counter(annotations[:, idx2]).values())
            total_averages[idx2] += largest_matching_count / annotations_length
        elif idx2 == 3:
            adjusted_annotations = [element.split('_')[0] for element in annotations[:, 2]]
            largest_matching_count = max(Counter(adjusted_annotations).values())
            total_averages[idx2] += largest_matching_count / annotations_length
        else:
            adjusted_annotations = [[*element][0] for element in annotations[:, 2]]
            largest_matching_count = max(Counter(adjusted_annotations).values())
            total_averages[idx2] += largest_matching_count / annotations_length

print(f"\n\nMatching clustid %: {(total_averages[0] / chains_length) * 100}")
print(f"Matching clustname %: {(total_averages[1] / chains_length) * 100}")
print(f"Matching subclass %: {(total_averages[2] / chains_length) * 100}")
print(f"Matching first term %: {(total_averages[3] / chains_length) * 100}")
print(f"Matching first letter %: {(total_averages[4] / chains_length) * 100}")

print(f"\nTime taken by script: {time.time() - time1}")
