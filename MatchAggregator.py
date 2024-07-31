import sys
import time
import numpy as np


def combine_matches(match_data_struct):
    # def build_chain(chain, index):
    #     built_chain = []
    #     for idx in range(index):
    #         built_chain.append(-99999)
    #     built_chain.append(chain[0])
    #     built_chain.append(chain[1])
    #     return built_chain
    #
    # all_chains = []
    # chain_map = {}  # Map to keep track of chains by their last element
    #
    # for idx, match_list in enumerate(match_data_struct):
    #     print(f"Chaining matches {idx}...")
    #     new_chains = []
    #     for a, b in match_list:
    #         if a in chain_map:
    #             chain_map[a] += (b,)
    #             chain_map[b] = chain_map[a]
    #             del chain_map[a]
    #         else:
    #             new_chain = build_chain((a, b), idx)
    #             new_chains.append(new_chain)
    #             chain_map[b] = new_chain
    #
    #     all_chains.extend(new_chains)
    #
    # for chain in chain_map:
    #     all_chains.append(chain)
    #
    # # Create chains without -99999
    # max_length = len(match_data_struct) + 1
    # valid_chains = [chain for chain in all_chains if isinstance(chain, (list, np.ndarray))
    #                 and -99999 not in chain
    #                 and len(chain) == max_length]
    #
    # return valid_chains

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

match_paths = [sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7]]
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
