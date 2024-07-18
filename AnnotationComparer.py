import numpy as np
import scanpy as sc
import sys

# Load AnnData objects
adata1_path = sys.argv[1]
adata2_path = sys.argv[2]
np_array_path = sys.argv[3]

np.set_printoptions(edgeitems=15)

adata1 = sc.read_h5ad(adata1_path)
adata2 = sc.read_h5ad(adata2_path)

print(adata1)
print(adata2)

matches = np.load(np_array_path)
print(matches)
print(f"Matches length: {len(matches)}")

annotations = ['clustid', 'clustname', 'subclass']

matched_clustid = 0
matched_clustname = 0
matched_subclass = 0

# Loop through the array and retrieve the cell ID and annotation
for idx, match in enumerate(matches):

    if sum(match) == 0: break

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
    print("\n")

print(f"Matching clustid %: {(matched_clustid / len(matches)) * 100}")
print(f"Matching clustname %: {(matched_clustname / len(matches)) * 100}")
print(f"Matching subclass %: {(matched_subclass / len(matches)) * 100}")