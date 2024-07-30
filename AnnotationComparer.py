import numpy as np
import scanpy as sc
import sys
import random


def permute_tuples(array):
    print("WARNING PERMUTE IS NOT WORKING!!!! THESE RESULTS ARE NOT VALID!!!!")
    # Extract the second elements of each tuple
    second_elements = [t[1] for t in array]

    # Shuffle the second elements randomly
    random.shuffle(second_elements)

    # Create new tuples with the original first elements and permuted second elements
    permuted_array = np.array([(array[i][0], second_elements[i]) for i in range(len(array))])

    return permuted_array


# Load AnnData objects
match_paths = [sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8]]
adata_struct = [sc.read(file) for file in match_paths]

match_chains = sys.argv[9]

np.set_printoptions(edgeitems=25)

print(adata_struct)

matches = np.load(match_chains)
# matches = np.array([t for t in matches if sum(t) != 0])
matches_length = len(matches)
print(f"Matches length: {matches_length}")

permute_matches = False
if permute_matches:
    matches = permute_tuples(matches)

print(matches)

annotations = ['clustid', 'clustname', 'subclass']

matched_clustid = 0
matched_clustname = 0
matched_subclass = 0
matched_subclass_first_term = 0
matched_first_letter = 0

# Loop through the array and retrieve the cell ID and annotation
for idx, match in enumerate(matches):

    # if sum(match) == 0: matches_length = idx + 1; break

    print(f"Match: {match}")

    annotations = adata_struct[:].obs_names[match[:]].obs[annotations].to_numpy()

    print(f"Annotations: {annotations}")
    if all(element == annotations[0, 0] for element in annotations[:, 0]): matched_clustid += 1
    if all(element == annotations[0, 1] for element in annotations[:, 1]): matched_clustname += 1
    if all(element == annotations[0, 2] for element in annotations[:, 2]): matched_subclass += 1
    if all(element.split('_')[0] == annotations[0, 2].split('_')[0] for element in annotations[:, 2]):
        matched_subclass_first_term += 1
    if all([*element][0] == [*(annotations[0, 2])][0] for element in annotations[:, 2]):
        matched_first_letter += 1

    print("\n")

print(f"Matching clustid %: {(matched_clustid / matches_length) * 100}")
print(f"Matching clustname %: {(matched_clustname / matches_length) * 100}")
print(f"Matching subclass %: {(matched_subclass / matches_length) * 100}")
print(f"Matching first term %: {(matched_subclass_first_term / matches_length) * 100}")
print(f"Matching first letter %: {(matched_first_letter / matches_length) * 100}")
