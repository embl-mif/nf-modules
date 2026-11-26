#!/usr/bin/env python3

"""
Load an intensity OME-Zarr and a label OME-Zarr with ZarrSmith, measure every
label with skimage.measure.regionprops (a base property set, plus any
user-requested structural and/or per-channel intensity properties), and
write the resulting feature table as a single CSV.

Per channel, the base intensity set is mean/min/max, the standard deviation
("std_intensity"), the coefficient of variation ("cv_intensity" = std/mean)
and the requested percentiles ("p25_intensity", "p75_intensity" by default)
of each label's pixel intensities.

Each label also gets a local background: the unlabelled (label 0) pixels
within --background-radius times its own equivalent radius of its edge.
Every label's ring is measured on its own, so the rings of nearby labels can
overlap and share background pixels; pixels of any label are never part of
a ring. Per channel, the
background median and standard deviation are written, together with the
signal-to-background ratio ("sbr_intensity" = mean / background median) and
the signal-to-noise ratio ("snr_intensity" = (mean - background median) /
background std).
"""

import argparse
import math
import os
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import zarrsmith.microscope_manager as mm
from scipy import ndimage
from skimage.measure import regionprops_table

# Shape/geometry properties always measured (when supported for the label
# image's dimensionality). "perimeter"/"orientation" are 2D-only in skimage
# and are used to derive "form_factor" (4*pi*area/perimeter^2).
BASE_STRUCTURAL_PROPERTIES_2D = (
    "label", "area", "bbox", "centroid", "perimeter", "orientation",
    "eccentricity", "solidity", "extent",
    "major_axis_length", "minor_axis_length",
)
BASE_STRUCTURAL_PROPERTIES_ND = (
    "label", "area", "bbox", "centroid", "solidity", "extent",
)

BASE_INTENSITY_PROPERTIES = ("label", "mean_intensity", "min_intensity", "max_intensity")

# Per-channel statistics always measured, keyed by every name (old and new
# skimage spelling) a user might also pass in --extra-intensity-properties,
# mapped to the column that already holds it. Requesting one of these again
# would duplicate a column, so it is skipped. The std is not part of
# BASE_INTENSITY_PROPERTIES: it is always added by measure_intensity through
# the pixel_std_intensity extra function.
ALREADY_MEASURED_INTENSITY_PROPERTIES = {
    "mean_intensity": "mean_intensity",
    "intensity_mean": "mean_intensity",
    "min_intensity": "min_intensity",
    "intensity_min": "min_intensity",
    "max_intensity": "max_intensity",
    "intensity_max": "max_intensity",
    "std_intensity": "std_intensity",
    "intensity_std": "std_intensity",
}

# Percentiles of each label's pixel intensities measured per channel when
# --intensity-percentiles is not given.
DEFAULT_INTENSITY_PERCENTILES = "25,75"

# Width, in multiples of each label's own equivalent radius, of the ring of
# unlabelled pixels around it used as its local background.
DEFAULT_BACKGROUND_RADIUS = 5.0


def pixel_std_intensity(regionmask: np.ndarray, intensity_image: np.ndarray) -> float:
    """Population standard deviation (ddof=0) of the label's pixel intensities.

    Named so it cannot collide with a skimage property: skimage >= 0.26 has a
    deprecated "std_intensity" alias that would shadow a custom function of
    that name. The column is renamed to "std_intensity" in measure_intensity.
    """
    return float(np.std(intensity_image[regionmask]))


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Measure skimage regionprops for every label of a label OME-Zarr, "
            "against every channel of an intensity OME-Zarr, and write the "
            "features to a CSV."
        )
    )
    parser.add_argument("--intensity", required=True, type=Path)
    parser.add_argument("--labels", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--intensity-resolution-level", type=int, default=0)
    parser.add_argument("--labels-resolution-level", type=int, default=0)
    parser.add_argument(
        "--channels",
        type=str,
        default=None,
        help=(
            "Comma/semicolon/whitespace-separated intensity channel indices "
            "to measure. Defaults to every channel."
        ),
    )
    parser.add_argument(
        "--extra-properties",
        type=str,
        default=None,
        help=(
            "Comma/semicolon/whitespace-separated extra skimage.measure."
            "regionprops properties to add to the base structural set "
            "(measured once per label, not per channel)."
        ),
    )
    parser.add_argument(
        "--extra-intensity-properties",
        type=str,
        default=None,
        help=(
            "Comma/semicolon/whitespace-separated extra skimage.measure."
            "regionprops properties that require an intensity image, "
            "measured once per requested channel (e.g. 'weighted_centroid')."
        ),
    )
    parser.add_argument(
        "--intensity-percentiles",
        type=str,
        default=DEFAULT_INTENSITY_PERCENTILES,
        help=(
            "Comma/semicolon/whitespace-separated percentiles (0-100) of each "
            "label's pixel intensities to measure per channel, written as "
            "'p<percentile>_intensity_<channel>' (e.g. '5,25,50,75,95'). "
            f"Defaults to '{DEFAULT_INTENSITY_PERCENTILES}'."
        ),
    )
    parser.add_argument(
        "--background-radius",
        type=float,
        default=DEFAULT_BACKGROUND_RADIUS,
        help=(
            "Width of the ring of unlabelled (label 0) pixels around each "
            "label used as its local background for the bg_*, sbr_intensity "
            "and snr_intensity columns, in multiples of the label's "
            "equivalent radius (the radius of a circle/sphere of the same "
            "area/volume). 0 disables these columns. Defaults to "
            f"{DEFAULT_BACKGROUND_RADIUS:g} cell radii."
        ),
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Threads used to measure the background rings. Defaults to the CPUs available to the process.",
    )
    parser.add_argument("--dataset-id", type=str, default=None)
    return parser.parse_args()


def available_cpus() -> int:
    # sched_getaffinity respects cgroup/cpuset limits (e.g. SLURM); it is not
    # available on every platform.
    if hasattr(os, "sched_getaffinity"):
        return len(os.sched_getaffinity(0))
    return os.cpu_count() or 1


def parse_token_list(text: str | None) -> list[str]:
    if not text:
        return []
    return [token for token in re.split(r"[;,\s]+", text.strip()) if token]


def parse_int_list(text: str | None) -> list[int]:
    tokens = parse_token_list(text)
    try:
        return [int(token) for token in tokens]
    except ValueError as exc:
        raise ValueError(
            f"Invalid channel list '{text}'. Expected integer channel "
            "indices separated by comma, semicolon, or whitespace."
        ) from exc


def parse_percentile_list(text: str | None) -> list[float]:
    tokens = parse_token_list(text)
    try:
        values = [float(token) for token in tokens]
    except ValueError as exc:
        raise ValueError(
            f"Invalid percentile list '{text}'. Expected numbers between 0 "
            "and 100 separated by comma, semicolon, or whitespace."
        ) from exc
    for value in values:
        if not 0 <= value <= 100:
            raise ValueError(f"Percentile {value:g} is out of range; expected a value between 0 and 100")
    return sorted(set(values))


def drop_already_measured(extra_intensity_properties: list[str]) -> list[str]:
    kept = []
    for prop in extra_intensity_properties:
        existing = ALREADY_MEASURED_INTENSITY_PROPERTIES.get(prop)
        if existing:
            print(f"Skipping extra intensity property '{prop}': already measured as {existing}")
        elif prop not in kept:
            kept.append(prop)
    return kept


def percentile_column_name(percentile: float) -> str:
    # 25 -> "p25_intensity", 2.5 -> "p2_5_intensity" (no "." so the name stays
    # usable in zarrsmith_filter's "property:min:max" expressions).
    text = np.format_float_positional(percentile, trim="-")
    return f"p{text.replace('.', '_')}_intensity"


def safe_divide(numerator: pd.Series, denominator: pd.Series) -> pd.Series:
    """numerator / denominator, NaN where the denominator is 0 or NaN."""
    return numerator / denominator.where(denominator != 0)


def equivalent_radius(area: np.ndarray, ndim: int) -> np.ndarray:
    """Radius of the ndim-ball whose volume is `area` pixels."""
    unit_ball_volume = math.pi ** (ndim / 2) / math.gamma(ndim / 2 + 1)
    return (area / unit_ball_volume) ** (1 / ndim)


def compact_labels(label_image: np.ndarray) -> np.ndarray:
    """Return the label image in an integer dtype of at most 4 bytes if possible.

    Integer labels of up to 4 bytes are kept as they are (no copy). Wider or
    non-integer labels are cast to int32 when their ids fit, else to int64.
    """
    if label_image.dtype.kind in "iu" and label_image.dtype.itemsize <= 4:
        return label_image
    max_label = int(label_image.max()) if label_image.size else 0
    target = np.int32 if max_label <= np.iinfo(np.int32).max else np.int64
    return label_image.astype(target, copy=False)


def measure_local_background(
    label_image: np.ndarray,
    channel_images: list[np.ndarray],
    regions: pd.DataFrame,
    radius_factor: float,
    threads: int = 1,
) -> tuple[pd.DataFrame, list[pd.DataFrame]]:
    """Measure each label's local background ring in every channel.

    A label's ring is every unlabelled (label 0) pixel within `radius_factor`
    times the label's equivalent radius of the label's edge. Each label is
    measured on its own, so the rings of nearby labels can overlap and share
    pixels. Pixels of any label are never part of a ring.

    `regions` is the structural regionprops table: its "label", "area" and
    "bbox-*" columns give each label's size and bounding box, so the label
    image is not scanned again. Labels are measured on `threads` threads
    (scipy's distance transform releases the GIL).

    Returns a table of ring sizes ("label", "bg_area") and, for each channel
    in `channel_images`, a table of "label", "bg_median_intensity" and
    "bg_std_intensity" (population std, ddof=0). Labels with an empty ring
    have bg_area 0 and NaN background intensities.
    """
    ndim = label_image.ndim
    labels = regions["label"].to_numpy()
    max_distances = radius_factor * equivalent_radius(regions["area"].to_numpy(dtype=np.float64), ndim)
    bbox_min = regions[[f"bbox-{axis}" for axis in range(ndim)]].to_numpy()
    bbox_max = regions[[f"bbox-{axis + ndim}" for axis in range(ndim)]].to_numpy()

    def measure_one(i: int) -> tuple[int, list[float], list[float]]:
        label, max_distance = labels[i], max_distances[i]
        margin = int(np.ceil(max_distance))
        crop = tuple(
            slice(max(int(lo) - margin, 0), min(int(hi) + margin, size))
            for lo, hi, size in zip(bbox_min[i], bbox_max[i], label_image.shape)
        )
        crop_labels = label_image[crop]
        # Distance of every pixel in the crop to the nearest pixel of this label.
        distance = ndimage.distance_transform_edt(crop_labels != label)
        ring = (crop_labels == 0) & (distance <= max_distance)
        del distance
        n_pixels = int(np.count_nonzero(ring))
        if not n_pixels:
            return 0, [np.nan] * len(channel_images), [np.nan] * len(channel_images)
        medians, stds = [], []
        for channel_image in channel_images:
            # Native dtype: np.median and np.std accumulate in float64 anyway.
            values = channel_image[crop][ring]
            medians.append(float(np.median(values)))
            stds.append(float(np.std(values, dtype=np.float64)))
        return n_pixels, medians, stds

    if threads > 1 and len(labels) > 1:
        with ThreadPoolExecutor(max_workers=threads) as executor:
            results = list(executor.map(measure_one, range(len(labels))))
    else:
        results = [measure_one(i) for i in range(len(labels))]

    area_df = pd.DataFrame({"label": labels, "bg_area": np.array([r[0] for r in results], dtype=np.int64)})
    channel_dfs = [
        pd.DataFrame(
            {
                "label": labels,
                "bg_median_intensity": [r[1][k] for r in results],
                "bg_std_intensity": [r[2][k] for r in results],
            }
        )
        for k in range(len(channel_images))
    ]
    return area_df, channel_dfs


def sanitize_suffix(text: str) -> str:
    return re.sub(r"[^0-9A-Za-z]+", "_", text).strip("_") or "channel"


def measure_structural_regions(label_image: np.ndarray, extra_properties: list[str]) -> pd.DataFrame:
    is_2d = label_image.ndim == 2
    base_properties = BASE_STRUCTURAL_PROPERTIES_2D if is_2d else BASE_STRUCTURAL_PROPERTIES_ND

    try:
        table = regionprops_table(label_image, properties=base_properties)
    except NotImplementedError:
        base_properties = BASE_STRUCTURAL_PROPERTIES_ND
        table = regionprops_table(label_image, properties=base_properties)

    df = pd.DataFrame(table)

    if "perimeter" in df.columns:
        with np.errstate(divide="ignore", invalid="ignore"):
            df["form_factor"] = np.where(
                df["perimeter"] > 0,
                4 * np.pi * df["area"] / np.square(df["perimeter"]),
                0.0,
            )

    if extra_properties:
        # Requested explicitly by the user: let a NotImplementedError (e.g. a
        # 2D-only property against a 3D label image) surface as a clear error
        # rather than being silently dropped like the opportunistic base set.
        extra_table = regionprops_table(label_image, properties=("label", *extra_properties))
        df = df.merge(pd.DataFrame(extra_table), on="label")

    return df


def measure_intensity(
    label_image: np.ndarray,
    intensity_image: np.ndarray,
    suffix: str,
    extra_intensity_properties: list[str],
    percentiles: list[float],
    background: pd.DataFrame | None,
) -> pd.DataFrame:
    properties = BASE_INTENSITY_PROPERTIES + tuple(extra_intensity_properties)
    extra_functions = [pixel_std_intensity]

    if percentiles:
        q = np.asarray(percentiles, dtype=float)

        def intensity_percentiles(regionmask, intensity_image):
            # One call per label computes every requested percentile at once.
            return np.percentile(intensity_image[regionmask], q)

        extra_functions.append(intensity_percentiles)

    table = regionprops_table(
        label_image,
        intensity_image=intensity_image,
        properties=properties,
        extra_properties=extra_functions,
    )
    df = pd.DataFrame(table)

    # regionprops_table expands an array-valued property into
    # "intensity_percentiles-0", "intensity_percentiles-1", ...
    df = df.rename(
        columns={
            "pixel_std_intensity": "std_intensity",
            **{
                f"intensity_percentiles-{i}": percentile_column_name(percentile)
                for i, percentile in enumerate(percentiles)
            },
        }
    )

    df["cv_intensity"] = safe_divide(df["std_intensity"], df["mean_intensity"])

    if background is not None:
        # A label with no unlabelled pixel in its ring (e.g. packed in by its
        # neighbours) has NaN background columns.
        df = df.merge(background, on="label", how="left")
        df["sbr_intensity"] = safe_divide(df["mean_intensity"], df["bg_median_intensity"])
        df["snr_intensity"] = safe_divide(df["mean_intensity"] - df["bg_median_intensity"], df["bg_std_intensity"])

    return df.rename(columns={col: f"{col}_{suffix}" for col in df.columns if col != "label"})


def extract_frame(level_data, indices: dict[int, int | None]) -> np.ndarray:
    indexer = [slice(None)] * level_data.ndim
    for axis, value in indices.items():
        if value is not None:
            indexer[axis] = value
    return np.asarray(level_data[tuple(indexer)])


def match_shape(array: np.ndarray, target_shape: tuple[int, ...]) -> np.ndarray:
    """Crop/zero-pad `array` at the high end of each axis to `target_shape`.

    Independently-built multiscale pyramids (e.g. ZarrSmith's own converter vs.
    ngff-zarr's to_multiscales, used by Cellpose) can round a coarse level's
    dimensions a pixel or two differently even from the same base image.
    Mirrors _match_mask_shape in cellpose_sam_zarr.py.
    """
    if array.shape == target_shape:
        return array
    if array.ndim != len(target_shape):
        raise ValueError(f"Rank mismatch: array shape {array.shape} vs target {target_shape}")

    adjusted = array
    for axis, target in enumerate(target_shape):
        current = adjusted.shape[axis]
        if current > target:
            slicer = [slice(None)] * adjusted.ndim
            slicer[axis] = slice(0, target)
            adjusted = adjusted[tuple(slicer)]
        elif current < target:
            pad_width = [(0, 0)] * adjusted.ndim
            pad_width[axis] = (0, target - current)
            adjusted = np.pad(adjusted, pad_width=pad_width, mode="constant", constant_values=0)

    return adjusted


def main():
    args = parse_args()

    channel_indices = parse_int_list(args.channels)
    extra_properties = parse_token_list(args.extra_properties)
    extra_intensity_properties = drop_already_measured(parse_token_list(args.extra_intensity_properties))
    percentiles = parse_percentile_list(args.intensity_percentiles)
    if args.background_radius < 0:
        raise ValueError(f"--background-radius must be >= 0, got {args.background_radius:g}")
    threads = args.threads if args.threads is not None else available_cpus()
    if threads < 1:
        raise ValueError(f"--threads must be >= 1, got {threads}")

    intensity_source = mm.ZarrManager(str(args.intensity))
    intensity_axes = intensity_source.metadata["axes"].lower()
    if args.intensity_resolution_level < 0 or args.intensity_resolution_level >= len(intensity_source.data):
        raise ValueError(
            f"--intensity-resolution-level {args.intensity_resolution_level} is out of "
            f"range for input with {len(intensity_source.data)} level(s)"
        )
    intensity_data = intensity_source.data[args.intensity_resolution_level]

    if "c" not in intensity_axes:
        raise ValueError(f"Intensity dataset does not contain a channel axis: axes={intensity_axes}")
    c_index = intensity_axes.index("c")
    n_channels = intensity_data.shape[c_index]

    if not channel_indices:
        channel_indices = list(range(n_channels))
    for channel_idx in channel_indices:
        if channel_idx < 0 or channel_idx >= n_channels:
            raise ValueError(f"Channel index {channel_idx} is out of range for {n_channels} channel(s)")

    channel_names = intensity_source.metadata.get("channel_names")
    channel_suffixes = {}
    for channel_idx in channel_indices:
        if channel_names and channel_idx < len(channel_names):
            channel_suffixes[channel_idx] = f"ch{channel_idx}_{sanitize_suffix(str(channel_names[channel_idx]))}"
        else:
            channel_suffixes[channel_idx] = f"ch{channel_idx}"

    labels_source = mm.ZarrManager(str(args.labels))
    labels_axes = labels_source.metadata["axes"].lower()
    if args.labels_resolution_level < 0 or args.labels_resolution_level >= len(labels_source.data):
        raise ValueError(
            f"--labels-resolution-level {args.labels_resolution_level} is out of "
            f"range for input with {len(labels_source.data)} level(s)"
        )
    labels_data = labels_source.data[args.labels_resolution_level]

    t_index_int = intensity_axes.index("t") if "t" in intensity_axes else None
    t_index_lbl = labels_axes.index("t") if "t" in labels_axes else None
    # A singleton T axis on the intensity image is dropped when labels have none.
    squeeze_int_t = (
        t_index_int is not None
        and t_index_lbl is None
        and intensity_data.shape[t_index_int] == 1
    )
    if squeeze_int_t:
        print(f"Dropping singleton T axis from intensity image (axes={intensity_axes})")
    elif (t_index_int is not None) != (t_index_lbl is not None):
        raise ValueError(
            "Intensity and labels disagree on whether a time axis is present "
            f"(intensity axes={intensity_axes}, labels axes={labels_axes})"
        )
    if (
        not squeeze_int_t
        and t_index_int is not None
        and intensity_data.shape[t_index_int] != labels_data.shape[t_index_lbl]
    ):
        raise ValueError(
            "Intensity and labels have a different number of timepoints "
            f"({intensity_data.shape[t_index_int]} vs {labels_data.shape[t_index_lbl]})"
        )

    if squeeze_int_t:
        timepoints = [None]
        int_t_fixed = {t_index_int: 0}
    else:
        timepoints = range(intensity_data.shape[t_index_int]) if t_index_int is not None else [None]
        int_t_fixed = {}

    frame_tables = []

    for t in timepoints:
        label_image = extract_frame(labels_data, {t_index_lbl: t} if t_index_lbl is not None else {})
        label_image = compact_labels(label_image)

        first_channel_image = extract_frame(
            intensity_data,
            {c_index: channel_indices[0], **int_t_fixed, **({t_index_int: t} if t is not None and t_index_int is not None else {})},
        )
        if label_image.shape != first_channel_image.shape:
            print(
                f"Adjusting label image shape from {label_image.shape} to "
                f"{first_channel_image.shape} (pyramid rounding mismatch between "
                f"--intensity-resolution-level={args.intensity_resolution_level} and "
                f"--labels-resolution-level={args.labels_resolution_level})"
            )
            label_image = match_shape(label_image, first_channel_image.shape)

        df = measure_structural_regions(label_image, extra_properties)

        def load_channel(channel_idx):
            if channel_idx == channel_indices[0]:
                return first_channel_image
            return extract_frame(
                intensity_data,
                {c_index: channel_idx, **int_t_fixed, **({t_index_int: t} if t is not None and t_index_int is not None else {})},
            )

        backgrounds = [None] * len(channel_indices)
        if args.background_radius > 0:
            # Every channel is held in memory so that each label's background
            # ring is computed once and measured in all channels.
            channel_images = [load_channel(channel_idx) for channel_idx in channel_indices]
            bg_area, backgrounds = measure_local_background(
                label_image, channel_images, df, args.background_radius, threads
            )
            df = df.merge(bg_area, on="label", how="left")
        else:
            # Without background columns, one channel at a time is enough.
            channel_images = None

        for k, (channel_idx, background) in enumerate(zip(channel_indices, backgrounds)):
            channel_image = channel_images[k] if channel_images is not None else load_channel(channel_idx)
            intensity_df = measure_intensity(
                label_image,
                channel_image,
                channel_suffixes[channel_idx],
                extra_intensity_properties,
                percentiles,
                background,
            )
            df = df.merge(intensity_df, on="label")

        if t is not None:
            df.insert(0, "timepoint", t)
        frame_tables.append(df)
        # Free this frame's images before the next frame is read.
        del label_image, first_channel_image, channel_images, channel_image

    features = pd.concat(frame_tables, ignore_index=True) if frame_tables else pd.DataFrame()
    if args.dataset_id:
        features.insert(0, "dataset_id", args.dataset_id)

    features.to_csv(args.output, index=False)

    intensity_source.close()
    labels_source.close()

    print(f"Measured channel(s): {[(idx, channel_suffixes[idx]) for idx in channel_indices]}")
    print(f"Wrote {len(features)} row(s) across {len(timepoints)} timepoint(s) to {args.output}")


if __name__ == "__main__":
    main()
