# resample_organo_cyte

Same as [organo_cyte](../organo_cyte/README.md), with a resampling step after conversion:
converts microscopy datasets to OME-Zarr, changes their pixel size, segments each
image tile by tile with Cellpose-SAM, merges labels across tile borders into one
multiscale label OME-Zarr, and measures the labels against the image, writing one
features CSV per dataset.

Segmentation and features run on the resampled image, so the labels have its pixel
size and every pixel-based setting (tile sizes, `cellpose_diameter`, ROI in voxels)
refers to the resampled grid. A row without a pixel size is not resampled: it is
segmented and measured on the converted image, exactly as in `organo_cyte`.

## Steps

1. `ZARRSMITH_CONVERSION` – convert the input to OME-Zarr (`zarrsmith convert`). Skipped
   if the store already exists.
2. `ZARRSMITH_RESAMPLE` – write a copy of the converted OME-Zarr with the new pixel
   size (`zarrsmith resample`). Only for rows with a `pixel_size` (or with
   `--pixel_size` set). Skipped if the store already exists.
3. `TILED_CELLPOSE_SAM_ZARR` (subworkflow) – segment the image tile by tile and merge
   the result into the final label OME-Zarr (`labels_output`):
   1. `TILE_CONSTRUCTOR` – compute the tile grid (and optional ROI).
   2. `PREALLOCATE_ZARR` – pre-allocate the merged `uint32` label store.
   3. `CELLPOSE_SAM_ZARR_MIF` – segment every tile on a GPU. With `-resume`,
      tiles already segmented in a previous run are reused.
   4. `FIND_TILE_NEIGHBORS` / `LOCAL_LABEL_MATCHER` – match labels between
      neighbouring tiles.
   5. `LABEL_SOLVER` – solve the label-matching graph into a global mapping.
   6. `WRITE_RELABELED_TILE` – write every tile into the merged label store.
   7. `POPULATE_PYRAMID` – fill pyramid levels 1..N from level 0 (nearest-neighbour
      `Sample`, so label ids are preserved).
4. `ZARRSMITH_FEATURE_EXTRACTION` – measure the labels against the (resampled) image
   (`skimage.measure.regionprops`) and write a features CSV.

## Input CSV

One dataset per row. Blank cells fall back to the default listed below. Only
`input` is mandatory. Paths containing commas must be CSV-quoted. Every output
column (`output`, `resampled_output`, `labels_output`, `features_output`)
must be an **absolute** path, or blank to use the default.

### General

| CSV column | Required | Default if not set | Description |
| --- | --- | --- | --- |
| `input` | **yes** | – | Path to the input microscopy dataset. |
| `dataset_id` | no | `<input basename>_<row number>` | Identifier used in task names and default output names. Keep it unique within the CSV. |

### 1. Conversion (`zarrsmith convert`)

Same columns as [organo_cyte](../organo_cyte/README.md#1-conversion-zarrsmith-convert):
`output`, `microscope`, `chunk_size`, `max_size(MB)`, `scene_index`, `zarr_format`,
`num_levels`, `downscale_factor`, `subset`, `channel_names`, `channel_colors`,
`shards`, `shard_target_mb`, `shard_exclude_axes`, `num_workers`, `drop_singleton`.
`output` defaults to `<dataset_id>.ome.zarr` in the work dir.

### 2. Resampling (`zarrsmith resample`)

| CSV column | Default if not set | Description |
| --- | --- | --- |
| `pixel_size` | `--pixel_size`, else no resampling | New pixel size in the dataset units (usually µm): one value for every spatial axis (`1`), or three in Z Y X order (`1 0.5 0.5`). |
| `resampled_output` | `<dataset_id>_resampled.ome.zarr` in the work dir | Path of the resampled OME-Zarr. Absolute path to keep it. Skipped if the store already exists. |
| `fast` | `false` | `true` to read from the coarsest pyramid level still fine enough instead of level 0: much less data when shrinking, but keeps the aliasing of that level. |
| `include_groups` | `true` | `false` to skip image subgroups. |
| `include_labels` | `true` | `false` to skip labels stored in the image. |

`num_workers` (writer threads) is shared with the conversion step.

The copy keeps the compression and number of levels of the converted image. Images
use linear interpolation with anti-aliasing on shrinking axes; labels use nearest
neighbour. See [ZARRSMITH_RESAMPLE](../../modules/local/zarrsmith_resample/README.md).

### 3. Tiled Cellpose-SAM segmentation

Same columns as [organo_cyte](../organo_cyte/README.md#2-tiled-cellpose-sam-segmentation):
`tile_size_x`, `tile_size_y`, `tile_size_z`, `tile_overlap`, `resolution_level`,
`use_physical_units`, `x_min` ... `z_max`, `cellpose_channels`, `cellpose_diameter`,
`cellpose_niter`, `cellpose_flow_threshold`, `cellpose_cellprob_threshold`,
`cellpose_halo`, `do_3D`, and `labels_output` (default
`<outdir>/masks/<dataset_id>_merged_mask.ome.zarr`).

### 4. Feature extraction

Same columns as [organo_cyte](../organo_cyte/README.md#3-feature-extraction):
`features_output` (default `<outdir>/measurements/<dataset_id>_features.csv`),
`intensity_resolution_level`, `labels_resolution_level`, `channels`,
`extra_properties`, `extra_intensity_properties`.

### Example

```csv
dataset_id,input,microscope,output,pixel_size,resampled_output,tile_size_x,tile_size_y,tile_overlap,cellpose_diameter,cellpose_channels,labels_output,channels,features_output
sample_01,/g/mif-hd/shared/data/sample_01.ome.tiff,opera,/scratch/user/zarr/sample_01.ome.zarr,1 0.5 0.5,/scratch/user/zarr/sample_01_iso.ome.zarr,2048,2048,128,30,0,/scratch/user/masks/sample_01_mask.ome.zarr,0;1,/scratch/user/features/sample_01.csv
sample_02,/g/mif-hd/shared/data/sample_02.ome.tiff,opera,,1,,,,,25,0,,,
```

## Parameters

| Parameter | Default | Purpose |
| --- | --- | --- |
| `--input_csv` | – (required) | Input CSV. |
| `--pixel_size` | none | Pixel size for rows with an empty `pixel_size`, e.g. `1` or `"1 0.5 0.5"`. Leave unset to resample only the rows that give one. |
| `--outdir` | `results` | Base directory for the default labels (`<outdir>/masks`) and features CSV (`<outdir>/measurements`). |
| `--resolution_level` | `0` | Pyramid level for tiling/segmentation when a row has none. |
| `--intensity_resolution_level` | `0` | Pyramid level for feature extraction on the image when a row has none. |
| `--labels_resolution_level` | `0` | Pyramid level for feature extraction on the mask when a row has none. |
| `--label_offset_step` | `100000` | Label-ID offset between tiles. Must exceed the largest per-tile label count. |
| `--multicut_beta` | `0.5` | `LABEL_SOLVER` edge-cut trade-off. |
| `--cellpose_model_path` | none | Custom pretrained Cellpose model. Defaults to Cellpose-SAM. |

## Running the workflow

```bash
module load Nextflow/25.10.1

nextflow run /path/to/nf-modules/workflows/resample_organo_cyte/main.nf \
    -lib /path/to/nf-modules/lib \
    --input_csv /path/to/input.csv \
    --cellpose_model_path /g/mif-hd/shared/cellpose_models/cpsam \
    -resume
```

For batch submission, update the marked values in `template_run.sh`, create
the log directory, and submit:

```bash
mkdir -p logs
sbatch template_run.sh
```

Per-tile segmentations are kept only in the Nextflow work directory; relaunch with
`-resume` to reuse finished tiles, and clean the work directory once the results
are checked.
