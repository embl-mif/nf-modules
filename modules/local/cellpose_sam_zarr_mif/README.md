# CELLPOSE_SAM_ZARR_MIF

Segments an OME-Zarr image, or one region of it (an ROI or a tile, plus a halo), with Cellpose-SAM (GPU) and writes the instance labels as an OME-Zarr label image. Retries on GPU out-of-memory / init failures (up to 2 times).

Reading and writing go through ZarrSmith (image `cellpose_zarr:4.1.1-zs1`, see `containers/cellpose_zarr`):

- **Input:** NGFF 0.4, 0.5 or 0.6, at `resolution_level`. The region is read with its halo (`ZarrManager.read_region`), Cellpose segments the block with the halo, and only the region itself is kept. The halo is only used with an ROI or tile bounds.
- **Output:** uint32 labels marked as label data, with the voxel size and units of the segmented level and the translation of the region's first voxel, so the store sits where the region is in the image. NGFF 0.4 inputs give 0.4 / zarr v2, 0.5 and 0.6 inputs give 0.5 / zarr v3.
- **Pyramid:** untiled outputs follow the source image from `resolution_level` on (same number of levels below it, same factors), sampled with `nearest` so ids are never mixed. Tiles (`meta.is_tiled`) are written at level 0 only, since they are only read there before they are fused.
- **Physical units:** with `use_physical_units`, ROI bounds and halo are converted by ZarrSmith with the NGFF 0.6 convention (the translation is the centre of the first voxel); a bound past the image edge is refused.

## Input
| Name | Type | Description |
| --- | --- | --- |
| `meta` | val | Sample map (`id`, plus Cellpose options such as `segmentation_channels`, `diameter`, `halo`, `resolution_level`, the ROI `x_min` ... `z_max`, `is_tiled`) |
| `local_image` | path | Local OME-Zarr to segment |
| `image_uri` | val | Path given as a value instead of `local_image` (local paths only: remote stores are not supported) |
| `pretrained_cellpose_model` | path | Optional custom Cellpose model file |

## Output
| Name | Type | Description |
| --- | --- | --- |
| `labels` | `tuple(meta, zarr)` | Label OME-Zarr (`meta.zarr_output`, default `<id>_labels.ome.zarr`) |
