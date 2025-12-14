"""Concatenate multiple AnnData objects into a single AnnData object."""

from pathlib import Path

import anndata as ad

# Define the directory path
adata_path = Path(
    "/Volumes/sep22/home/wet_lab/_Experiments/009_ST_Xenium/out/data/adata"
)

# Get all .h5ad files
files = sorted(adata_path.glob("*.h5ad"))
print(f"Found {len(files)} files to concatenate.")

# Load all adata objects
adatas = []
for file in files:
    adata = ad.read_h5ad(file)
    # Add sample ID - this is crucial for later visualization
    adata.obs["sample_id"] = file.stem
    adatas.append(adata)

# Concatenate all objects
adata_combined = ad.concat(
    adatas,
    join="outer",  # Keep all genes
    label="batch",  # Creates 'batch' column
    keys=[f.stem for f in files],
    index_unique="_",
)

# Save the combined object
output_path = adata_path / "combined_adata.h5ad"

adata_combined.write_h5ad(output_path)
print(f"Combined AnnData saved to {output_path}")
