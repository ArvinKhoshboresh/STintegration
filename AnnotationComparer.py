import numpy as np
import scanpy as sc
import sys
import random

def jaccard_similarity(str1, str2):
    # Remove all "_" characters from the strings
    str1 = str1.replace("_", "")
    str2 = str2.replace("_", "")

    # Convert the strings to sets of characters
    set1 = set(str1)
    set2 = set(str2)

    # Calculate the intersection and union of the sets
    intersection = set1.intersection(set2)
    union = set1.union(set2)

    # Calculate the Jaccard similarity
    similarity = len(intersection) / len(union)

    return similarity


def permute_tuples(array):
    # Extract the second elements of each tuple
    second_elements = [t[1] for t in array]

    # Shuffle the second elements randomly
    random.shuffle(second_elements)

    # Create new tuples with the original first elements and permuted second elements
    permuted_array = np.array([(array[i][0], second_elements[i]) for i in range(len(array))])

    return permuted_array


# Load AnnData objects
adata1_path = sys.argv[1]
adata2_path = sys.argv[2]
np_array_path = sys.argv[3]

np.set_printoptions(edgeitems=25)

adata1 = sc.read_h5ad(adata1_path)
adata2 = sc.read_h5ad(adata2_path)

print(adata1)
print(adata2)

matches = np.load(np_array_path)
matches = np.array([t for t in matches if sum(t) != 0])
matches_length = len(matches)
print(f"Matches length: {matches_length}")
permute_matches = False
if permute_matches: matches = permute_tuples(matches)
print(matches)

annotations = ['clustid', 'clustname', 'subclass']

matched_clustid = 0
matched_clustname = 0
matched_subclass = 0
matched_subclass_first_term = 0
matched_first_letter = 0
jaccard_tracker = 0

# Loop through the array and retrieve the cell ID and annotation
for idx, match in enumerate(matches):

    if sum(match) == 0: matches_length = idx + 1; break

    print(f"Idx: {idx}    Match: {match}")
    cell_index1 = match[0]
    cell_index2 = match[1]

    annotation1 = adata1[str(cell_index1)].obs[annotations].to_numpy()[0]
    annotation2 = adata2[str(cell_index2)].obs[annotations].to_numpy()[0]

    print(f'{cell_index1}: Annotations: {annotation1}')
    print(f'{cell_index2}, Annotations: {annotation2}')
    if annotation1[0] == annotation2[0]: matched_clustid += 1
    if annotation1[1] == annotation2[1]: matched_clustname += 1
    if annotation1[2] == annotation2[2]: matched_subclass += 1
    if annotation1[2].split('_')[0] == annotation2[2].split('_')[0]: matched_subclass_first_term += 1
    if [*annotation1[2]][0] == [*annotation2[2]][0]: matched_first_letter += 1
    jaccard_tracker += jaccard_similarity(annotation1[2], annotation2[2])

    print("\n")

print(f"Matching clustid %: {(matched_clustid / matches_length) * 100}")
print(f"Matching clustname %: {(matched_clustname / matches_length) * 100}")
print(f"Matching subclass %: {(matched_subclass / matches_length) * 100}")
print(f"Matching first term %: {(matched_subclass_first_term / matches_length) * 100}")
print(f"Matching first letter %: {(matched_first_letter / matches_length) * 100}")
print(f"Average Jaccard index: {(jaccard_tracker / matches_length) * 100}")