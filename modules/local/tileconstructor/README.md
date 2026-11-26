# TILE_CONSTRUCTOR

Splits an OME-Zarr image into a regular grid of (optionally overlapping) tiles and writes their bounds, in voxels as half-open `[min, max)` intervals, to a CSV. Each row has a unique `dataset_id` of the form `<source_dataset_id>__<tile_id>`. Tiles are numbered from 1, x fastest, then y, then z.

The tiles are planned by ZarrSmith (`zarrsmith.utils.tiling.plan_tiles`, branch `consolidated_translation`, not released yet), from the image as ZarrSmith reads it, so they agree with ZARRSMITH_FUSE_LABELS.

- An axis given `tile_size_<ax>` is cut at that size.
- An axis without one is cut only when a tile would read more than `max_tile_voxels` voxels (halo included, times `tile_resize ** ndim`). ZarrSmith then chooses the number of tiles along each such axis that fits with the least cut area, then the fewest tiles. If every axis has a tile size and the tiles are over the budget, the task fails.
- Without `tile_size_<ax>` and without `max_tile_voxels`, the axis is one tile (the store's chunk shape is no longer used).

## Input
| Name | Type | Description |
| --- | --- | --- |
| `meta` | val | Sample map, see below |
| `local_image` | path | Source OME-Zarr (or a `NO_*` sentinel when given by URI) |
| `image_uri` | val | URI of the source OME-Zarr, or empty |

### `meta` keys
| Key | Description |
| --- | --- |
| `id` | Source image name (`source_dataset_id`) |
| `tile_size_x`, `tile_size_y`, `tile_size_z` | Tile size per axis. Default: not split, unless `max_tile_voxels` cuts it |
| `max_tile_voxels` | Most voxels one tile may read, halo included |
| `tile_resize` | Linear factor the next step enlarges tiles by (e.g. `30 / diameter` for Cellpose). Default: 1 |
| `tile_overlap` | Overlap between neighbouring tiles. Default: 0 |
| `halo` | Context the next step reads around each tile; counts towards `max_tile_voxels` and is the smallest last tile before it is merged into the previous one |
| `resolution_level` | Pyramid level to tile. Default: 0 |
| `x_min` ... `z_max` | Region of interest. Default: the whole level |
| `use_physical_units` | Sizes, overlap, halo and ROI in physical units (the NGFF translation is the centre of the first voxel). Default: voxels |

## Output
| Name | Type | Description |
| --- | --- | --- |
| `tiles` | `tuple(meta, csv)` | Tile table (`*_tiles.csv`). `tile_size_<ax>`/`tile_overlap_<ax>` are the planned values; an axis left whole reports its ROI length and overlap 0 |

Each row also names its forward neighbours from the tile plan: `right_dataset_id`/`right_tile_index` (next tile along x), `bottom_*` (along y) and `behind_*` (along z), blank at the end of an axis. LOCAL_LABEL_MATCHER compares each pair of touching tiles once, from the lower tile's side.
