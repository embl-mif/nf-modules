# ORGANOIDCYTE_FEATURE_EXTRACTION

Measures every label of a label OME-Zarr against an intensity OME-Zarr with `skimage.measure.regionprops` (shape properties plus per-channel intensities) and writes one CSV row per label.

It is an extended version of [`ZARRSMITH_FEATURE_EXTRACTION`](../zarrsmith_feature_extraction/README.md): it measures the same shape and intensity properties, and also measures the spread of each label's pixel intensities in every channel (see [Additional features](#additional-features)).

## Input
| Name | Type | Description |
| --- | --- | --- |
| `meta` | val | Sample map (`id`, `features_*` options and output name) |
| `intensity_zarr` | path | Intensity OME-Zarr |
| `labels_zarr` | path | Label OME-Zarr |

## Output
| Name | Type | Description |
| --- | --- | --- |
| `csv` | `tuple(meta, csv)` | Feature table (`meta.csv_output`, default `<id>_features.csv`) |
| `versions` | path | `versions.yml` |

## Measured features

Shape properties are measured once per label:

- 2D labels: `area`, `bbox`, `centroid`, `perimeter`, `orientation`, `eccentricity`, `solidity`, `extent`, `major_axis_length`, `minor_axis_length`, and the derived `form_factor` (`4*pi*area/perimeter^2`).
- 3D labels: `area`, `bbox`, `centroid`, `solidity`, `extent`.
- Any `meta.features_extra_properties`.

Intensity properties are measured once per channel in `meta.features_channels` (every channel by default):

- `mean_intensity`, `min_intensity`, `max_intensity` (the same as `ZARRSMITH_FEATURE_EXTRACTION`).
- `std_intensity` and the percentile columns (new, see below).
- Any `meta.features_extra_intensity_properties`.

Each intensity column gets a per-channel suffix: `<property>_ch<index>`, or `<property>_ch<index>_<channel name>` when the OME-Zarr has channel names. For example, `std_intensity_ch0` or `p75_intensity_ch1_GFP`.

## Additional features

`ZARRSMITH_FEATURE_EXTRACTION` does not compute these columns. They are added for every measured channel:

| Column | Description |
| --- | --- |
| `std_intensity_<channel>` | Population standard deviation (`ddof=0`) of the label's pixel intensities. Always measured. |
| `p<percentile>_intensity_<channel>` | Percentiles of the label's pixel intensities, computed with NumPy's default linear interpolation. By default these are the 25th and 75th (`p25_intensity`, `p75_intensity`). |

Percentiles are set with `meta.features_intensity_percentiles` (the `--intensity-percentiles` option of the script). Give values from 0 to 100, separated by `;`, `,` or whitespace, for example `5;25;50;75;95`. Duplicates are removed and the values are sorted. A decimal point in the column name becomes `_`: `2.5` gives `p2_5_intensity`. This keeps the name usable in `ZARRSMITH_FILTER` `property:min:max` expressions. A value outside 0–100 raises an error.

If `meta.features_extra_intensity_properties` asks for a property that is already measured, it is skipped so that the column is not duplicated. This covers `mean_intensity`, `min_intensity`, `max_intensity` and `std_intensity`, in either skimage spelling (`intensity_mean`, etc.).

## `meta` options

| Key | Default if unset | Description |
| --- | --- | --- |
| `id` | – | Dataset id, written to the `dataset_id` column. |
| `csv_output` | `<id>_features.csv` | Output CSV path. If it already exists, the extraction is skipped. |
| `features_intensity_resolution_level` | `0` | Pyramid level read from the intensity image. |
| `features_labels_resolution_level` | `0` | Pyramid level read from the labels. |
| `features_channels` | every channel | Intensity channel indices to measure. |
| `features_extra_properties` | none | Extra shape `regionprops` properties. |
| `features_extra_intensity_properties` | none | Extra intensity `regionprops` properties (e.g. `weighted_centroid`). |
| `features_intensity_percentiles` | `25;75` | Percentiles to measure per channel. |
