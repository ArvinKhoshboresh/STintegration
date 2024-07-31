import sys
import random
import time

import numpy as np


def combine_matches(match_data_struct):
    def extend_chain(chain, length):
        while len(chain) < length:
            chain += (-99999,)
        return chain

    def build_chain(chain, index):
        built_chain = ()
        for idx in range(index):
            built_chain += (-99999,)
        built_chain += chain
        return built_chain

    all_chains = []

    for idx, match_list in enumerate(match_data_struct):
        print(f"Chaining matches {idx}...")
        new_chains = []
        for a, b in match_list:
            found_chain = False
            for chain in all_chains:
                if chain[-1] == a:
                    chain += b
                    found_chain = True
                    break
            if not found_chain:
                new_chain = build_chain((a, b), idx)
                new_chains.append(new_chain)

        all_chains.extend(new_chains)

    # Ensure all chains are filled with -99999 up to the maximum length
    max_length = len(match_data_struct) + 1
    all_chains = [extend_chain(chain, max_length) for chain in all_chains]

    # Create chains without -99999
    valid_chains = [chain for chain in all_chains if -99999 not in chain]

    return valid_chains, all_chains




    def extend_chain(chain, length):
        return chain + (-99999,) * (length - len(chain))

    def build_chain(chain, index):
        return (-99999,) * index + chain

    all_chains = []
    chain_map = {}  # Map to keep track of chains by their last element

    for idx, match_list in enumerate(match_data_struct):
        print(f"Chaining matches {idx}...")
        new_chains = []
        for a, b in match_list:
            if a in chain_map:
                chain_map[a] += (b,)
                chain_map[b] = chain_map[a]
                del chain_map[a]
            else:
                new_chain = build_chain((a, b), idx)
                new_chains.append(new_chain)
                chain_map[b] = new_chain

        all_chains.extend(new_chains)

    # Ensure all chains are filled with -99999 up to the maximum length
    max_length = len(match_data_struct) + 1
    all_chains = [extend_chain(chain, max_length) for chain in all_chains]

    # Create chains without -99999
    valid_chains = [chain for chain in all_chains if -99999 not in chain]

    return valid_chains, all_chains


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


def remove_trailing_zeros(tuples_list):
    reversed_list = tuples_list[::-1]

    while (reversed_list[0] == (0, 0)).all():
        reversed_list.pop(0)

    return reversed_list[::-1]


time1 = time.time()

np.set_printoptions(edgeitems=10)

permute_matches = False
if permute_matches:
    match_path = sys.argv[1]
    matches = np.load(match_path)
    match_data_struct = [permute_all_tuples(matches), permute_all_tuples(matches), permute_all_tuples(matches),
                         permute_all_tuples(matches), permute_all_tuples(matches), permute_all_tuples(matches),
                         permute_all_tuples(matches), permute_all_tuples(matches)]
else:
    match_paths = [sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7]]
    match_data_struct = [np.load(file) for file in match_paths]

match_data_struct = remove_trailing_zeros(match_data_struct[:])

# print(match_files)

valid_chains, all_chains = combine_matches(match_data_struct)

print(all_chains)

print(f"Total number of full chains: {len(valid_chains)}\n")

print(valid_chains)

# Save valid matches to disk if needed
np.save('all_chains.npy', np.array(all_chains))
np.save('valid_chains.npy', np.array(valid_chains))

print(f"Script Took {time.time() - time1}s")
