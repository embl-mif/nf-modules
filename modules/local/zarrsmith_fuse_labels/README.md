# ZARRSMITH_FUSE_LABELS

Merges the label tiles of one tiled segmentation (CELLPOSE_SAM_ZARR_MIF, one store per tile) into one OME-Zarr label image, with `zarrsmith fuse-labels --tiles_csv ... --label_offset_step ... --label_mapping ...` (no blending, label data). It replaces PREALLOCATE_ZARR, WRITE_RELABELED_TILE and POPULATE_PYRAMID in the tiled Cellpose subworkflow.

- **Placement:** each tile sits at the NGFF translation of its own store, so no tile bounds are needed.
- **Numbering:** labels are numbered as LOCAL_LABEL_MATCHER numbered the graph nodes (every non-zero id + `tile_index * label_offset_step`, from the tiles CSV; zarrsmith first checks that the staged stores and the CSV list the same tiles, and writes the offsets it used to `<output name>.label_offsets.csv` next to the output), then LABEL_SOLVER's global mapping (`--label_mapping`) merges the objects split between tiles.
- **No blending:** every voxel is copied from one tile, and an overlap is cut in its middle. Pyramid levels sample voxels instead of averaging them, so ids are never mixed.
- **Pyramid:** follows the source image from `resolution_level` on (`--infer_downscale_factors <source> --start_level <resolution_level>`; `num_levels`/`downscale_factor` override it), so level k of the labels matches level `resolution_level + k` of the image.
- **Layout:** chunks, shards and the number of fusion workers follow `zarrsmith fuse`'s defaults (about 4 MB chunks, 256 MB shards, workers cut to what fits in memory, at most `task.cpus`).

The output is written in full each time (an existing store is replaced, once the fusion has succeeded). zarrsmith also writes `<output>.transforms.json` (the position of every tile) next to the output.

Requires ZarrSmith with `zarrsmith fuse-labels` (branch `consolidated_translation`, not released yet).

## Input
| Name | Type | Description |
| --- | --- | --- |
| `meta` | val | Sample map, see below |
| `tile_labels` | path | Label stores of every tile, staged in `tiles/` |
| `tiles_csv` | path | `<id>_tiles_and_segmentations.csv`: `tile_index` and `segmented_tile_path` of every tile |
| `global_mapping` | path | LABEL_SOLVER's global label mapping (JSON) |
| `local_image` | path | Source image (or a `NO_*` sentinel when given by URI), for the pyramid levels |
| `image_uri` | val | Source image URI, or empty |

### `meta` keys
| Key | Required | Description |
| --- | --- | --- |
| `id` | yes | Source image name |
| `label_offset_step` | yes | Same value as LOCAL_LABEL_MATCHER |
| `zarr_output` | no | Output path. Default: `<id>_merged_labels.ome.zarr` in the work directory |
| `resolution_level` | no | Source level that was segmented. Default: 0 |
| `num_levels` | no | Pyramid levels. Default: as the source from `resolution_level` on |
| `downscale_factor` | no | Factor between levels, one value or `z y x`. Default: as the source |
| `dtype` | no | Output dtype. Default: `uint32` |
| `pyramid_method` | no | `nearest` or `mode` for the pyramid levels. Default: `nearest` (zarrsmith's default for labels) |

Extra options of `zarrsmith fuse-labels` (e.g. `-cs`, `-ms`, `--shards`, `-cp`) can be passed with `ext.args`.

## Output
| Name | Type | Description |
| --- | --- | --- |
| `zarr` | tuple(val(meta), path) | Merged OME-Zarr label image |
| `versions` | path | `versions.yml` |
