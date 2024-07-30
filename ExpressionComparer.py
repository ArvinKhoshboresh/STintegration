import numpy as np
import pandas as pd
from statsmodels.multivariate.manova import MANOVA
import sys
import scanpy as sc
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

# Loop through the array and retrieve the cell ID and annotation
for idx, match in enumerate(matches):
    print(f"Match: {match}")

    X_vectors = adata_struct[:].obs_names[match[:]].X.to_numpy()

    print(X_vectors)

# Convert lists of vectors into a single DataFrame
data = np.vstack([list1, list2])
labels = ['group1'] * len(list1) + ['group2'] * len(list2)

# Create a DataFrame with the data and labels
df = pd.DataFrame(data, columns=['var1', 'var2', 'var3'])
df['group'] = labels

# Perform MANOVA
maov = MANOVA.from_formula('var1 + var2 + var3 ~ group', data=df)
print(maov.mv_test())
