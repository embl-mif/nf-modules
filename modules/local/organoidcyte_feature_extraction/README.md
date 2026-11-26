# ORGANOIDCYTE_FEATURE_EXTRACTION

Measures every label of a label OME-Zarr against an intensity OME-Zarr with `skimage.measure.regionprops` (shape properties plus per-channel intensities) and writes one CSV row per label.

It is an extended version of [`ZARRSMITH_FEATURE_EXTRACTION`](../zarrsmith_feature_extraction/README.md): it measures the same shape and intensity properties. It also measures, in every channel, the spread of each label's pixel intensities and how bright the label is compared with its local background (see [Additional features](#additional-features)).

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
- `bg_area`: the number of pixels in the label's background ring (see [Local background](#local-background)). Not written when `meta.features_background_radius` is `0`.
- Any `meta.features_extra_properties`.

Intensity properties are measured once per channel in `meta.features_channels` (every channel by default):

- `mean_intensity`, `min_intensity`, `max_intensity` (the same as `ZARRSMITH_FEATURE_EXTRACTION`).
- `std_intensity`, `cv_intensity`, the percentile columns and the background columns (new, see below).
- Any `meta.features_extra_intensity_properties`.

Each intensity column gets a per-channel suffix: `<property>_ch<index>`, or `<property>_ch<index>_<channel name>` when the OME-Zarr has channel names. For example, `std_intensity_ch0` or `p75_intensity_ch1_GFP`.

## Additional features

`ZARRSMITH_FEATURE_EXTRACTION` does not compute these columns. They are added for every measured channel:

| Column | Description |
| --- | --- |
| `std_intensity_<channel>` | Population standard deviation (`ddof=0`) of the label's pixel intensities. Always measured. |
| `cv_intensity_<channel>` | Coefficient of variation: `std_intensity / mean_intensity`. Always measured. |
| `p<percentile>_intensity_<channel>` | Percentiles of the label's pixel intensities, computed with NumPy's default linear interpolation. By default these are the 25th and 75th (`p25_intensity`, `p75_intensity`). |
| `bg_median_intensity_<channel>` | Median of the label's background ring pixels. |
| `bg_std_intensity_<channel>` | Population standard deviation (`ddof=0`) of the label's background ring pixels. This is the noise used for `snr_intensity`. |
| `sbr_intensity_<channel>` | Signal-to-background ratio: `mean_intensity / bg_median_intensity`. |
| `snr_intensity_<channel>` | Signal-to-noise ratio: `(mean_intensity - bg_median_intensity) / bg_std_intensity`. |

Percentiles are set with `meta.features_intensity_percentiles` (the `--intensity-percentiles` option of the script). Give values from 0 to 100, separated by `;`, `,` or whitespace, for example `5;25;50;75;95`. Duplicates are removed and the values are sorted. A decimal point in the column name becomes `_`: `2.5` gives `p2_5_intensity`. This keeps the name usable in `ZARRSMITH_FILTER` `property:min:max` expressions. A value outside 0–100 raises an error.

A ratio whose denominator is `0` is written as empty (NaN), for example `cv_intensity` of a label whose mean is `0`.

### Local background

A label's background is the ring of unlabelled (label `0`) pixels around it, up to `meta.features_background_radius` times the label's own equivalent radius away from its edge (default `5`, the `--background-radius` option of the script). The equivalent radius is the radius of a circle (2D) or sphere (3D) with the same area or volume as the label, so a larger cell gets a wider ring. For example, a cell 30 pixels across (radius 15) gets a ring 75 pixels wide with the default `5`. Distances are in pixels at the labels' resolution level, the same along every axis, so on anisotropic 3D data the ring is thicker in z than in xy in physical units.

Each label's ring is measured on its own: a background pixel within reach of several labels is part of all their rings, even if it is closer to another label. Pixels that belong to any label, the label itself or its neighbours, are never part of a ring.

Each label's ring is measured separately, so this step takes longer with more labels and wider rings, especially in 3D. The labels are split over the task's CPUs (`task.cpus`, passed as `--threads`), so giving the process more CPUs makes it faster. While the rings are measured, all measured channels of a frame are held in memory together. With `meta.features_background_radius` set to `0`, the channels are read one at a time.

The median is used for the background level because it is less affected by bright unsegmented objects in the ring than the mean. The noise is the standard deviation of the ring, which is where there is no cell signal.

A label that is surrounded by other labels on every side has no ring: `bg_area` is `0` and its background, `sbr_intensity` and `snr_intensity` columns are empty. Use `bg_area` to drop labels whose ring has too few pixels for reliable values. Set the radius to `0` to skip all the background columns.

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
| `features_background_radius` | `5` | Width of each label's local background ring, in multiples of the label's equivalent radius. `0` disables the background columns. |
