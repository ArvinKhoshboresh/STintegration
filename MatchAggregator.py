import sys
import random
import numpy as np


def read_match_file(filename):
    return np.load(filename)


def combine_matches(match_data):

    combined_matches = {}
    for idx, match_array in enumerate(match_data):
        for match in match_array:
            start, end = match
            if start not in combined_matches:
                combined_matches[start] = [-999999] * 8
            combined_matches[start][idx] = end

    valid_matches = []
    for key, values in combined_matches.items():
        if all(value != -999999 for value in values):
            valid_matches.append((key, *values))

    return valid_matches


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


np.set_printoptions(edgeitems=100)

# match_files = [sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],
#                sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8]]
# match_data = [read_match_file(file) for file in match_files]

permute_matches = True
if permute_matches:
    np_array_path = sys.argv[1]
    matches = np.load(np_array_path)
    match_files = [permute_all_tuples(matches), permute_all_tuples(matches), permute_all_tuples(matches),
                   permute_all_tuples(matches), permute_all_tuples(matches), permute_all_tuples(matches),
                   permute_all_tuples(matches), permute_all_tuples(matches)]

# print(match_files)

valid_matches = combine_matches(match_files)
print("Total number of tuples with 8 valid matches:", len(valid_matches))

print(valid_matches)

# Save valid matches to disk if needed
# np.save('valid_matches.npy', np.array(valid_matches))
