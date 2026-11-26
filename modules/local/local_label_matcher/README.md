# LOCAL_LABEL_MATCHER

Builds one tile's contribution to the region-adjacency graph used for multicut stitching: boundary edges between its labels and the labels of its right / bottom / behind neighbours (label ids are offset per tile so they are globally unique).

## Input
| Name | Type | Description |
| --- | --- | --- |
| `meta` | val | Sample map (`id`, tile index, overlaps) |
| `self_zarr` | path | This tile's label zarr |
| `right_zarr`, `bottom_zarr`, `behind_zarr` | path | Neighbour label zarrs (placeholders when the neighbour does not exist) |

Each label zarr is read at level 0, the first dataset of its NGFF `multiscales` (NGFF 0.4 or 0.5): stores written by ngff-zarr (`scale0/nuclei_labels`) and by ZarrSmith (`0`) both work, and a store without that metadata is read at `scale0/nuclei_labels`. Non-spatial axes (`t`, `c`) are dropped by name and must have size 1. Only the overlap band of each neighbour is read.

## Output
| Name | Type | Description |
| --- | --- | --- |
| `local_graph` | `tuple(meta, json)` | Local graph (`*.json`), input of `LABEL_SOLVER` |
