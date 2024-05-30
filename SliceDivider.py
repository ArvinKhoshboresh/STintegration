import numpy as np
import anndata
import os

# Load AnnData object
adata = anndata.read_h5ad("C:/Users/arkho/OneDrive/Desktop/adata3.h5ad")

# Check if 'X_spatial' is in adata.obsm and extract z coordinates
if 'X_spatial' in adata.obsm:
    spatial_coords = adata.obsm['X_spatial']
    # Assuming spatial coordinates are stored in the order [x, y, z]
    z_coords = spatial_coords[:, 2]
else:
    raise KeyError("The key 'X_spatial' is not found in adata.obsm.")

# Define range for z coordinates and the step size
z_min, z_max = np.min(z_coords), np.max(z_coords)
step_size = 0.1

output_dir = 'z_splits'
os.makedirs(output_dir, exist_ok=True)

# Iterate through z coordinate ranges and create new AnnData objects
z_start = z_min
while z_start < z_max:
    z_end = z_start + step_size

    # Subset AnnData object
    mask = (z_coords >= z_start) & (z_coords < z_end)
    adata_subset = adata[mask, :].copy()

    # Print index and coordinates of each cell in the subset
    print(f'Subset z_{z_start:.1f}_to_{z_end:.1f}:')
    for idx in adata_subset.obs.index:
        coord = adata_subset.obsm['X_spatial'][adata_subset.obs.index.get_loc(idx)]
        print(f'  Cell: {idx}, Coordinates: {coord}')

    # Write subset to disk
    subset_filename = f'z_{z_start:.1f}_to_{z_end:.1f}.h5ad'
    subset_filepath = os.path.join(output_dir, subset_filename)
    adata_subset.write_h5ad(subset_filepath)

    # Update the start of the next range
    z_start = z_end

print("Splitting completed and files saved.")
