"""Convert Xenium output to AnnData format and save."""

from pathlib import Path

from spatialdata_io import xenium

sample_path = Path(
    "/Volumes/sep22/home/wet_lab/_Experiments/009_ST_Xenium/data/xenium/run_1/20251001__141239__SP25164_SARA_PATTI_RUN_1/output-XETG00431__0021047__COPD_R035_V3__20251001__141533"
)
roi = "COPD_R035_V3"
run_num = 1

sdata = xenium(sample_path)

adata_path = Path(
    "/Volumes/sep22/home/wet_lab/_Experiments/009_ST_Xenium/out/data/adata"
)
adata_path.mkdir(parents=True, exist_ok=True)

adata = sdata.tables["table"]
adata.obs["ROI"] = roi
adata.obs["run"] = run_num
adata.write_h5ad(adata_path / f"{roi}.h5ad")

print(f"Finished {roi}")
