import numpy as np
import scanpy as sc
import sys
import pandas as pd

# Load AnnData objects
adata1_path = sys.argv[1]
adata2_path = sys.argv[2]
np_array_path = sys.argv[3]
annotations_path = sys.argv[4]

adata1 = sc.read_h5ad(adata1_path)
adata2 = sc.read_h5ad(adata2_path)

print(adata1)
print(adata2)

matches = np.load(np_array_path)

annotations_df = pd.read_csv(annotations_path, index_col=0)

print(annotations_df)

annotations = ['neurotransmitter', 'class', 'subclass', 'supertype', 'cluster']

# Loop through the array and retrieve the cell ID and annotation
for idx, value in enumerate(matches):
    try:
        cell_id1 = adata1.obs_names[idx]  # Cell ID at the matches index position
        cell_id2 = adata2.obs_names[value]  # Cell ID at the matches value position

        row1 = annotations_df.loc[cell_id1]
        row2 = annotations_df.loc[cell_id2]

        annotation1 = row1[annotations]
        annotation2 = row2[annotations]

        print(f'Match: {idx}, {value}')
        print(f'{idx}: {cell_id1}, Annotations: {annotation1.to_dict()}')
        print(f'{value}: {cell_id2}, Annotations: {annotation2.to_dict()}')
        print("\n")
    except KeyError as e:
        print(f'Error: {e}')
