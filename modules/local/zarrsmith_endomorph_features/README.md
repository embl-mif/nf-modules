# ZARRSMITH_ENDOMORPH_FEATURES

Measures per-cell morphology (area, perimeter, aspect ratio, orientation) and the mean intensity of every measured channel in each band R1..RN and over the whole cell, the centroid and area of every gap not touching the image border (label `N+1` of the regions zarr; `N+3` gaps are excluded), and writes the per-image summary (cell and gap counts and empty-area statistics, read from the regions zarr). Makes no assumption about what the channels stain. Skipped if all three CSVs already exist.

## Input
| Name | Type | Description |
| --- | --- | --- |
| `meta` | val | Sample map (`id`, `endomorph_ring_depths` (only the number of bands is used), `endomorph_measure_channels` (default: all channels), resolution levels, output names) |
| `intensity_zarr` | path | Intensity OME-Zarr (the MIP) |
| `labels_zarr` | path | Final cell labels |
| `regions_zarr` | path | Bands + empty-area labels from `ZARRSMITH_ENDOMORPH_REGIONS` |

## Output
| Name | Type | Description |
| --- | --- | --- |
| `csv` | `tuple(meta, cells_csv, summary_csv, gaps_csv)` | Per-cell CSV, per-image summary CSV and per-gap CSV (`meta.gaps_csv_output`, default `<id>_endomorph_gaps.csv`) |
| `versions` | path | `versions.yml` |
