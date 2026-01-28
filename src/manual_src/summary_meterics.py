"""Script to concatenate metrics_summary.csv files from multiple output directories."""

from pathlib import Path

import pandas as pd

# Define variables


# Define the base directory
base_dir = Path(
    "/Volumes/phenotypingsputumasthmaticsaurorawellcomea1/live/Sara_Patti/009_ST_Xenium/data/xenium_raw/20251028__162706__SP25164_SARA_PATTI_RUN_2"
)

output_dir = Path(
    "/Volumes/sep22/home/wet_lab/_Experiments/009_ST_Xenium/data/out_data/summary_metrics/csv"
)

run = "_".join(base_dir.name.split("_")[-2:])

output_dir.mkdir(parents=True, exist_ok=True)

# List to store dataframes
dfs = []

# Iterate through directories starting with "output-"
#! Need to come back so it is dynamic and will go into run_2,
# !run_3 etc and look for the folder and files
for folder in base_dir.glob("output-*"):
    if folder.is_dir():
        csv_path = folder / "metrics_summary.csv"

        if csv_path.exists():
            print(f"Reading: {csv_path}")
            df = pd.read_csv(csv_path)

            # Add a column to track which folder the data came from
            df["source_folder"] = folder.name

            dfs.append(df)
        else:
            print(f"Warning: metrics_summary.csv not found in {folder.name}")

# Concatenate all dataframes
if dfs:
    # Concatenate by shared columns (intersection of all columns)
    combined_df = pd.concat(dfs, axis=0, ignore_index=True, join="inner")

    # Extract Condition from region_name
    if "region_name" in combined_df.columns:
        combined_df["Condition"] = combined_df["region_name"].str.split("_").str[0]
        print(
            f"\n Condition column added: {combined_df['Condition'].unique().tolist()}"
        )
    else:
        print("\nWarning: 'region_name' column not found")

    print(f"\nSuccessfully concatenated {len(dfs)} files")
    print(f"Combined dataframe shape: {combined_df.shape}")
    print(f"\nColumns in combined dataframe:\n{combined_df.columns.tolist()}")

    # Save the combined dataframe
    output_path = base_dir / f"combined_metrics_summary_{run}.csv"
    combined_df.to_csv(output_path, index=False)
    print(f"\nCombined data saved to: {output_path}")

    # Display first few rows
    print("\nFirst few rows:")
    print(combined_df.head())
else:
    print("No CSV files found to concatenate")

# Add column with run identifier
combined_df["run"] = int(str(base_dir).split("_")[-1])

# Save the combined dataframe
output_path = output_dir / f"combined_metrics_summary_{run}.csv"
combined_df.to_csv(output_path, index=False)
print(f"\nCombined data saved to: {output_path}")
