"""Subset the HLCA to a random selection of ~50k cells for testing purposes."""

import numpy as np
import scanpy as sc

input_path = "/rds/general/user/sep22/ephemeral/recode_hlca_full_processed.h5ad"
output_path = "/rds/general/user/sep22/ephemeral/recode_hlca_random_subset.h5ad"

# Load
adata = sc.read_h5ad(input_path)

# Pick ~50k cells (adjust to whatever your machine can handle)
target_n = 50000

if adata.n_obs > target_n:
    idx = np.random.choice(adata.n_obs, target_n, replace=False)
    adata_sub = adata[idx].copy()
else:
    adata_sub = adata.copy()

# Save
adata_sub.write(output_path)
print(f"Saved subset: {adata_sub.shape} → {output_path}")
