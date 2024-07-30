import sys
import random
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


def permute_all_tuples(array):
    # Extract the second elements of each tuple
    first_elements = [t[0] for t in array]
    second_elements = [t[1] for t in array]

    # Shuffle the second elements randomly
    random.shuffle(first_elements)
    random.shuffle(second_elements)

    # Create new tuples with the original first elements and permuted second elements
    permuted_array = np.array([(first_elements[i], second_elements[i]) for i in range(len(array))])

    return permuted_array


def remove_trailing_zeros(tuples_list):
    reversed_list = tuples_list[::-1]

    while reversed_list and reversed_list[0] == (0, 0):
        reversed_list.pop(0)

    return reversed_list[::-1]


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
# np.save('valid_matches.npy', np.array(valid_matches))
