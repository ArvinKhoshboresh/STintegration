import glob
import os
import sys
import time
import numpy as np


def combine_matches(match_data_struct):

    chains = match_data_struct[0]
    for idx, match_list in enumerate(match_data_struct[1:]):
        print(f"Chaining matches {idx + 1}...")
        chain_dict = {chain[-1]: chain for chain in chains}

        found_chains = []
        for a, b in match_list:
            try:
                new_chain = list(chain_dict[a])
                new_chain.append(b)
                found_chains.append(new_chain)
                print(new_chain)
            except:
                continue

        chains = found_chains

    return chains

def permute_all_tuples(matrix):
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


def remove_trailing_zeros(tuples_array):
    last_non_zero_index = len(tuples_array)
    for i in range(len(tuples_array) - 1, -1, -1):
        if not np.array_equal(tuples_array[i], [0, 0]):
            last_non_zero_index = i + 1
            break

    # Return the array up to the last non-(0, 0) index
    return tuples_array[:last_non_zero_index]


time1 = time.time()

# Path to the directory containing the .npy files
folder_path = 'Matches/'
match_paths = glob.glob(os.path.join(folder_path, '*.npy'))
match_data_struct = [np.load(file) for file in match_paths]
match_data_struct = [remove_trailing_zeros(matches) for matches in match_data_struct]

print(f"Num matches 0: {len(match_data_struct[0])}")

cut_data = False
if cut_data:
    np.random.seed(42)
    print("Data Split.")

    for idx, match_file in enumerate(match_data_struct):
        cut_data_factor = 1000
        num_chains = match_file.shape[0]
        indices = np.random.permutation(num_chains)
        split = indices[:num_chains // cut_data_factor]
        match_data_struct[idx] = match_file[split].copy()

    print(f"Num CUT matches 0: {len(match_data_struct[0])}")

valid_chains = combine_matches(match_data_struct)
# print(valid_chains)

print(f"Total number of full chains: {len(valid_chains)}\n")

# Save valid matches to disk if needed
np.save('valid_chains.npy', valid_chains)

print(f"Script Took {time.time() - time1}s")
