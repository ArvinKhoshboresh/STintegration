import numpy as np
import scanpy as sc
import logging
import sys
import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()

# Load AnnData objects
adata1_path = sys.argv[1]
adata2_path = sys.argv[2]
np_array_path = sys.argv[3]
annotations_path = sys.argv[4]

adata1 = sc.read_h5ad(adata1_path)
adata2 = sc.read_h5ad(adata2_path)

logger.info(adata1)

matches = np.load(np_array_path)

annotations_df = pd.read_csv(annotations_path, index_col=0)

# Loop through the array and retrieve the cell ID and annotation
for idx, value in enumerate(matches):
    try:
        cell_id1 = adata1.obs_names[idx]  # Cell ID at the index position
        cell_id2 = adata2.obs_names[value]  # Cell ID at the value position

        annotation1 = annotations_df.loc[
            cell_id1, ['neurotransmitter', 'class', 'subclass', 'supertype', 'cluster']]
        annotation2 = annotations_df.loc[
            cell_id2, ['neurotransmitter', 'class', 'subclass', 'supertype', 'cluster']]

        print(f'Match: {idx}, {value}')
        print(f'{idx}: {cell_id1}, Annotations: {annotation1.to_dict()}')
        print(f'{value}: {cell_id2}, Annotations: {annotation2.to_dict()}')
    except KeyError as e:
        print(f'Error: {e}')
