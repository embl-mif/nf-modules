#!/usr/bin/env python3

"""Segment an OME-Zarr image (or one region of it) with Cellpose-SAM and write the
labels as an OME-Zarr label image (CELLPOSE_SAM_ZARR_MIF).

Reading and writing go through ZarrSmith:

* the region (an ROI or one tile, plus its halo) is read with
  ``ZarrManager.read_region`` at the chosen resolution level, from NGFF 0.4, 0.5 or
  0.6 inputs;
* Cellpose segments the region with its halo, and only the region itself is kept;
* the labels are written with ``ArrayManager(..., data_type="label")``: uint32, the
  level's voxel size and units, and the translation of the region's first voxel, so
  the store sits where the region is in the image. NGFF 0.4 inputs give 0.4 / zarr v2
  outputs, 0.5 and 0.6 inputs give 0.5 / zarr v3 (``--ome-zarr-version`` overrides).

Pyramid: by default the labels follow the source image from the segmented level on
(same number of levels below it, same factors), sampled with "nearest" so ids are never
mixed. ``--num-levels 1`` writes level 0 only (tiles of a tiled segmentation, which are
only read at level 0 before they are fused).

Example:
    cellpose_sam_zarr.py -i image.ome.zarr -o labels.ome.zarr --channels 1 --diameter 30
"""

import argparse
import math
import os
import re
from typing import Any

import dask.array as da
import numpy as np
import zarr

from zarrsmith.microscope_manager import ArrayManager, ZarrManager
from zarrsmith.utils.dataset_ops import infer_downscale_factors

SPATIAL = ("z", "y", "x")


# ---------------------------------------------------------------------------
# Command-line parsing helpers
# ---------------------------------------------------------------------------

def _parse_index_cli_values(values: list[str] | None, option_name: str) -> list[int] | None:
    if values is None:
        return None
    tokens: list[str] = []
    for raw in values:
        text = str(raw).strip()
        if not text:
            continue
        text = re.sub(r'^"+|"+$', '', text)
        text = re.sub(r"^'+|'+$", '', text)
        tokens.extend(part for part in re.split(r"[;,\s]+", text) if part)
    if not tokens:
        return None
    try:
        return [int(token) for token in tokens]
    except ValueError as exc:
        raise ValueError(
            f"Invalid {option_name} value. Expected integer indices separated by semicolon, comma, or whitespace."
        ) from exc


def _parse_optional_float_cli(value: str) -> float | None:
    text = str(value).strip()
    return None if text.lower() == "none" else float(text)


def _parse_optional_int_cli(value: str) -> int | None:
    text = str(value).strip()
    return None if text.lower() == "none" else int(text)


def _parse_bool_cli(value: str) -> bool:
    text = str(value).strip().lower()
    if text in {"true", "1", "yes"}:
        return True
    if text in {"false", "0", "no", "none"}:
        return False
    raise argparse.ArgumentTypeError(f"invalid boolean value: {value!r}")


# ---------------------------------------------------------------------------
# Input: OME-Zarr version, timepoints, channels, region
# ---------------------------------------------------------------------------

def output_format(input_zarr: str, requested: str | None) -> tuple[str, int]:
    """(NGFF version, zarr format) of the output: as the input (0.6 is written as 0.5,
    which ZarrSmith writes), unless ``requested``."""
    if requested is not None:
        version = str(requested).strip()
    else:
        attrs = zarr.open_group(input_zarr, mode="r").attrs.asdict()
        ome = attrs.get("ome")
        if isinstance(ome, dict):
            version = "0.5"  # 0.5 and 0.6 keep their metadata under "ome"
        else:
            version = str((attrs.get("multiscales") or [{}])[0].get("version", "0.4"))
    if version not in {"0.4", "0.5"}:
        raise ValueError(f"Unsupported output OME-Zarr version '{version}'. Supported values are: 0.4, 0.5")
    return version, (2 if version == "0.4" else 3)


def _select_timepoints(shape: dict[str, int], timepoint_indices: list[int] | None) -> list[int | None]:
    if "t" not in shape:
        if timepoint_indices is not None:
            raise ValueError("Input has no time dimension, but --timepoints was provided")
        return [None]
    if timepoint_indices is None:
        return list(range(shape["t"]))
    selected = [int(t) for t in timepoint_indices]
    if not selected:
        raise ValueError("--timepoints was provided but no indices were specified")
    for t in selected:
        if not 0 <= t < shape["t"]:
            raise ValueError(f"timepoint index {t} out of range for {shape['t']} timepoint(s)")
    return selected


def _select_channels(shape: dict[str, int], channel_indices: list[int] | None) -> list[int] | None:
    if "c" not in shape:
        if channel_indices is not None:
            if channel_indices == [0] or channel_indices == 0:  # the only channel
                return None
            raise ValueError(f"Input has no channel dimension, but --channels was provided: {channel_indices}")
        return None
    if channel_indices is None:
        return list(range(shape["c"]))
    selected = [int(c) for c in channel_indices]
    if not selected:
        raise ValueError("--channels was provided but no channel indices were specified")
    for c in selected:
        if not 0 <= c < shape["c"]:
            raise ValueError(f"channel index {c} out of range for {shape['c']} channel(s)")
    return selected


def region_selection(
    axes: str,
    bounds: dict[str, tuple[float | None, float | None]],
    halo: float,
    use_physical_units: bool,
) -> tuple[dict[str, tuple], float | int]:
    """``read_region`` selection and halo from the --<ax>-min/--<ax>-max bounds.

    Voxel bounds are widened to whole voxels (floor/ceil); physical bounds are converted
    by ZarrSmith (every voxel overlapping them). The halo is only used with an ROI, as
    before: without bounds the whole level is read.
    """
    if "x" not in axes or "y" not in axes:
        raise ValueError(f"Input dims must contain x and y axes, got {list(axes)}")
    if halo < 0.0:
        raise ValueError(f"--halo must be >= 0, got {halo}")
    selection: dict[str, tuple] = {}
    for ax, (lo, hi) in bounds.items():
        if lo is None and hi is None:
            continue
        if ax not in axes:
            raise ValueError(f"Input has no {ax} dimension, but --{ax}-min/--{ax}-max was provided")
        if lo is not None and hi is not None and hi <= lo:
            raise ValueError(f"{ax}-bounds are invalid: min={lo} must be < max={hi}")
        if not use_physical_units:
            lo = None if lo is None else int(math.floor(lo + 1e-9))
            hi = None if hi is None else int(math.ceil(hi - 1e-9))
        selection[ax] = (lo, hi)
    if not selection:
        return {}, 0
    return selection, (halo if use_physical_units else int(math.ceil(halo - 1e-9)))


# ---------------------------------------------------------------------------
# Cellpose
# ---------------------------------------------------------------------------

def resolve_cellpose_gpu_device() -> dict[str, Any]:
    """Require a CUDA GPU and fail fast when GPU initialization is not possible."""
    from cellpose import core

    try:
        import torch
    except ImportError as exc:
        raise RuntimeError("GPU_INIT_ERROR: torch is not installed in this environment") from exc
    try:
        cuda_available = bool(torch.cuda.is_available())
        cuda_device_count = int(torch.cuda.device_count()) if cuda_available else 0
    except Exception as exc:
        raise RuntimeError(f"GPU_INIT_ERROR: failed to query CUDA availability: {exc}") from exc
    cuda_visible_devices = os.getenv("CUDA_VISIBLE_DEVICES", "<unset>")
    if not cuda_available or cuda_device_count < 1:
        raise RuntimeError(
            "GPU_INIT_ERROR: CUDA GPU not available "
            f"(is_available={cuda_available}, device_count={cuda_device_count}, "
            f"CUDA_VISIBLE_DEVICES={cuda_visible_devices})"
        )
    try:
        _ = torch.cuda.get_device_name(0)
    except Exception as exc:
        raise RuntimeError(f"GPU_INIT_ERROR: CUDA device discovery failed despite visible device(s): {exc}") from exc
    if not core.use_gpu():
        raise RuntimeError("GPU_INIT_ERROR: Cellpose reported no usable CUDA GPU despite CUDA availability")
    return {"gpu": True}


def load_cellpose_model(cellpose_model_path: str = ""):
    """The Cellpose model on the GPU (the pretrained Cellpose-SAM unless a model path is given)."""
    from cellpose import io, models

    io.logger_setup()
    model_kwargs = resolve_cellpose_gpu_device()
    print(f"Cellpose device settings: {model_kwargs}")
    if cellpose_model_path:
        return models.CellposeModel(pretrained_model=cellpose_model_path, **model_kwargs)
    return models.CellposeModel(**model_kwargs)


def _match_mask_shape(mask: np.ndarray, expected_shape: tuple[int, ...]) -> np.ndarray:
    """Crop or zero-pad ``mask`` to ``expected_shape``, at the high end of each axis
    (the origin stays where it is)."""
    if mask.ndim != len(expected_shape):
        raise ValueError(f"Cellpose mask rank mismatch: mask.ndim={mask.ndim}, expected ndim={len(expected_shape)}")
    adjusted = mask
    for axis, target in enumerate(expected_shape):
        current = adjusted.shape[axis]
        if current > target:
            slicer = [slice(None)] * adjusted.ndim
            slicer[axis] = slice(0, target)
            adjusted = adjusted[tuple(slicer)]
        elif current < target:
            pad_width = [(0, 0)] * adjusted.ndim
            pad_width[axis] = (0, target - current)
            adjusted = np.pad(adjusted, pad_width=pad_width, mode="constant", constant_values=0)
    return np.asarray(adjusted)


def segment_block(
    model,
    block: np.ndarray,
    block_axes: list[str],
    channels: list[int] | None,
    anisotropy: float,
    *,
    diameter: float | None,
    flow_threshold: float,
    cellprob_threshold: float,
    do_3D: bool,
    flow3D_smooth: float,
    stitch_threshold: float,
    niter: int | None,
    min_size: int,
) -> np.ndarray:
    """Cellpose labels of one timepoint of a block (axes ``block_axes`` without t),
    with the spatial shape of the block."""
    axes = list(block_axes)
    if "c" in axes:
        c_axis = axes.index("c")
        if channels is None:
            raise ValueError("No channel selection available despite channel dimension in input")
        if len(channels) == 1:
            block = np.take(block, channels[0], axis=c_axis)
            axes.pop(c_axis)
        else:
            block = np.take(block, channels, axis=c_axis)

    has_channel_axis = "c" in axes
    if has_channel_axis and block.ndim not in (3, 4):
        raise ValueError(f"Expected 2D/3D data with channel axis for Cellpose, got shape {block.shape}")
    if not has_channel_axis and block.ndim not in (2, 3):
        raise ValueError(f"Expected 2D or 3D data for Cellpose, got shape {block.shape}")
    n_slices = block.shape[axes.index("z")] if "z" in axes else 1

    eval_kwargs: dict[str, Any] = {
        "flow_threshold": flow_threshold,
        "cellprob_threshold": cellprob_threshold,
        "min_size": min_size,
        "niter": niter,
    }
    if diameter is not None:
        eval_kwargs["diameter"] = diameter
    if has_channel_axis:
        eval_kwargs["channel_axis"] = axes.index("c")
    if n_slices > 1:
        eval_kwargs["z_axis"] = axes.index("z")
        if do_3D:
            eval_kwargs.update(do_3D=True, anisotropy=anisotropy, flow3D_smooth=flow3D_smooth)
        else:
            eval_kwargs.update(do_3D=False, stitch_threshold=stitch_threshold)

    mask, _, _ = model.eval(block, **eval_kwargs)

    # Cellpose-SAM may return a mask a few voxels off the input shape.
    mask_axes = [ax for ax in axes if ax != "c"]
    expected = tuple(int(block.shape[axes.index(ax)]) for ax in mask_axes)
    mask = np.asarray(mask)
    if tuple(mask.shape) != expected:
        print(f"Adjusting Cellpose mask shape from {tuple(mask.shape)} to {expected}")
        mask = _match_mask_shape(mask, expected)
    return mask


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_labels(
    output_zarr: str,
    labels: np.ndarray,
    label_axes: str,
    md: dict,
    level: int,
    translation: dict[str, float],
    chunks: tuple[int, ...],
    ngff_version: str,
    zarr_format: int,
    num_levels: int | None,
) -> None:
    """Write ``labels`` as an OME-Zarr label image at ``translation``, with the voxel size
    and units of ``level`` of the source (``md``), and its pyramid."""
    spatial_in = [ax for ax in md["axes"] if ax in SPATIAL]
    spatial_out = [ax for ax in label_axes if ax in SPATIAL]
    level_scale = dict(zip(spatial_in, md["scales"][level]))
    units_in = md.get("units") or ()
    units = dict(zip(spatial_in, units_in)) if len(units_in) == len(spatial_in) else {}

    out_md = {
        "axes": label_axes,
        "size": [tuple(labels.shape)],
        "chunksize": [chunks],
        "dtype": "uint32",
        "data_type": "label",
        "scales": [tuple(float(level_scale[ax]) for ax in spatial_out)],
        "units": tuple(units.get(ax, "micrometer") for ax in spatial_out),
        "translations": [{ax: float(translation[ax]) for ax in spatial_out}],
    }
    if "t" in label_axes:
        for key in ("time_increment", "time_increment_unit"):
            if md.get(key) is not None:
                out_md[key] = md[key]

    levels, factor = infer_downscale_factors(md, level)  # the source's pyramid from `level` on
    if num_levels is not None:
        levels = max(1, int(num_levels))
    manager = ArrayManager([da.from_array(labels, chunks=chunks)], out_md)
    if levels > 1:
        manager.build_pyramid(num_levels=levels, downscale_factor=factor, method="nearest")
    print(f"Writing {levels} level(s) (factor {factor if levels > 1 else '-'}), NGFF {ngff_version}, "
          f"zarr v{zarr_format}, chunks {chunks}")
    # drop_singleton=False keeps the axes as segmented (e.g. a t of 1 or a z of 1).
    manager.to_zarr(output_zarr, ngff_version=ngff_version, zarr_format=zarr_format,
                    drop_singleton=False, overwrite=True, progress=False)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def cellpose_zarr(
    input_zarr: str,
    output_zarr: str,
    cellpose_model_path: str = "",
    channel_indices: list[int] | None = None,
    timepoint_indices: list[int] | None = None,
    x_min: float | None = None,
    x_max: float | None = None,
    y_min: float | None = None,
    y_max: float | None = None,
    z_min: float | None = None,
    z_max: float | None = None,
    halo: float = 0.0,
    use_physical_units: bool = False,
    diameter: float | None = None,
    flow_threshold: float = 0.4,
    cellprob_threshold: float = 0.0,
    do_3D: bool = False,
    flow3D_smooth: float = 0.0,
    stitch_threshold: float = 0.1,
    niter: int | None = None,
    min_size: int = 15,
    resolution_level: int = 0,
    ome_zarr_version: str | None = None,
    num_levels: int | None = None,
    model=None,
) -> None:
    print(f"Opening {input_zarr}")
    zm = ZarrManager(input_zarr, mode="r")
    md = zm.metadata
    n_levels = len(md["size"])
    if not 0 <= resolution_level < n_levels:
        raise ValueError(f"--resolution-level {resolution_level} is out of range for input with {n_levels} level(s)")
    axes = "".join(md["axes"])
    shape = dict(zip(axes, md["size"][resolution_level]))
    spatial = [ax for ax in axes if ax in SPATIAL]
    scale = dict(zip(spatial, md["scales"][resolution_level]))
    anisotropy = scale.get("z", 1.0) / scale["x"] if scale.get("x") else 1.0
    ngff_version, zarr_format = output_format(input_zarr, ome_zarr_version)

    timepoints = _select_timepoints(shape, timepoint_indices)
    channels = _select_channels(shape, channel_indices)
    selection, halo_arg = region_selection(
        axes, {"x": (x_min, x_max), "y": (y_min, y_max), "z": (z_min, z_max)}, halo, use_physical_units)

    if model is None:
        model = load_cellpose_model(cellpose_model_path)

    masks, region = [], None
    for t in timepoints:
        sel = dict(selection, **({"t": t} if t is not None else {}))
        region = zm.read_region(resolution_level, halo=halo_arg, physical=use_physical_units, **sel)
        block = region.data
        block_axes = list(region.axes)
        if "t" in block_axes:  # one timepoint at a time
            block = np.take(block, 0, axis=block_axes.index("t"))
            block_axes.remove("t")
        mask = segment_block(
            model, block, block_axes, channels, anisotropy,
            diameter=diameter, flow_threshold=flow_threshold, cellprob_threshold=cellprob_threshold,
            do_3D=do_3D, flow3D_smooth=flow3D_smooth, stitch_threshold=stitch_threshold,
            niter=niter, min_size=min_size,
        )
        # Segmented with the halo; keep the region itself.
        mask_axes = [ax for ax in block_axes if ax != "c"]
        crop = dict(zip(region.axes, region.crop))
        masks.append(mask[tuple(crop[ax] for ax in mask_axes)])

    label_axes = "".join((["t"] if "t" in axes else []) + mask_axes)
    labels = np.stack(masks) if "t" in axes else masks[0]
    if labels.size and labels.max() > np.iinfo(np.uint32).max:
        raise ValueError(f"{int(labels.max())} labels do not fit in uint32")
    labels = labels.astype(np.uint32)

    print(f"Processed resolution level: {resolution_level}")
    print(f"Processed timepoint(s): {timepoints}")
    print(f"Used channel(s): {channels}")
    print(f"Region (voxel): { {ax: region.bounds[ax] for ax in spatial} }")
    print(f"Region (physical): {region.physical_bounds('core')}")
    print(f"Read with halo (voxel): { {ax: region.window[ax] for ax in spatial} }")
    print(f"Halo read (voxel, before/after): {region.halo if selection else 'none (no ROI)'}")

    chunks_in = dict(zip(axes, md["chunksize"][resolution_level])) if md.get("chunksize") else {}
    chunks = tuple(min(int(chunks_in.get(ax, n)), int(n)) for ax, n in zip(label_axes, labels.shape))
    print(f"Writing labels to {output_zarr}")
    write_labels(output_zarr, labels, label_axes, md, resolution_level, region.core_translation(), chunks,
                 ngff_version, zarr_format, num_levels)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("-i", "--input_zarr", type=str, required=True)
    parser.add_argument("-o", "--output_zarr", type=str, required=True)
    parser.add_argument(
        "-c", "--channels", "--segmentation_channels", dest="segmentation_channels", type=str, nargs="+",
        default=None,
        help="Channel indices to process. Accepts semicolon/comma/space-separated values. If omitted, all channels are used.",
    )
    parser.add_argument("-m", "--cellpose_model_path", type=str, default="")
    parser.add_argument("-t", "--timepoints", type=int, nargs="+", default=None,
                        help="Timepoint indices to process. If omitted, all timepoints are used.")
    parser.add_argument("--ome-zarr-version", type=str, choices=["0.4", "0.5"], default=None,
                        help="OME-Zarr version of the output. Default: as the input (0.6 inputs give 0.5).")
    for ax in "xyz":
        parser.add_argument(f"--{ax}-min", type=float, default=None, help=f"Optional {ax}-axis ROI lower bound (inclusive).")
        parser.add_argument(f"--{ax}-max", type=float, default=None, help=f"Optional {ax}-axis ROI upper bound (exclusive).")
    parser.add_argument(
        "--halo", type=float, default=0.0,
        help="Optional non-negative ROI halo margin in the same units as --x-*/--y-*/--z-* (voxel by default, "
             "physical with --use-physical-units). Ignored when no ROI bounds are provided.",
    )
    parser.add_argument("--use-physical-units", action="store_true",
                        help="Interpret --x-*/--y-*/--z-* bounds and --halo in physical units instead of voxel units.")
    parser.add_argument("--resolution-level", type=int, default=0,
                        help="Zero-based OME-Zarr resolution level to operate on. Default: 0.")
    parser.add_argument("--num-levels", type=int, default=None,
                        help="Pyramid levels of the output. Default: as the source from --resolution-level on; "
                             "1 writes level 0 only (tiles of a tiled segmentation).")
    parser.add_argument("--diameter", type=_parse_optional_float_cli, default=None,
                        help="Expected cell diameter in pixels, or 'None' to auto-estimate.")
    parser.add_argument("--flow_threshold", type=float, default=0.4, help="Cellpose flow error threshold.")
    parser.add_argument("--cellprob_threshold", type=float, default=0.0, help="Cellpose cell probability threshold.")
    parser.add_argument("--do_3D", type=_parse_bool_cli, default=False,
                        help="Force 3D segmentation. If False (default), z-stacks are segmented plane by plane and stitched.")
    parser.add_argument("--flow3D_smooth", type=float, default=0.0, help="Sigma for smoothing flows in 3D segmentation.")
    parser.add_argument(
        "--stitch_threshold", type=float, default=0.1,
        help="IoU threshold used to stitch per-slice 2D masks into 3D labels when running non-3D segmentation on a z-stack.",
    )
    parser.add_argument("--niter", type=_parse_optional_int_cli, default=None,
                        help="Number of mask reconstruction iterations, or 'None' to use Cellpose's default.")
    parser.add_argument("--min_size", type=int, default=15, help="Minimum number of pixels for a valid mask.")
    return parser


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()
    try:
        segmentation_channels = _parse_index_cli_values(args.segmentation_channels, "--segmentation_channels/--channels")
    except ValueError as exc:
        parser.error(str(exc))
    cellpose_zarr(
        args.input_zarr, args.output_zarr, args.cellpose_model_path, segmentation_channels, args.timepoints,
        args.x_min, args.x_max, args.y_min, args.y_max, args.z_min, args.z_max, args.halo, args.use_physical_units,
        diameter=args.diameter, flow_threshold=args.flow_threshold, cellprob_threshold=args.cellprob_threshold,
        do_3D=args.do_3D, flow3D_smooth=args.flow3D_smooth, stitch_threshold=args.stitch_threshold,
        niter=args.niter, min_size=args.min_size, resolution_level=args.resolution_level,
        ome_zarr_version=args.ome_zarr_version, num_levels=args.num_levels,
    )


if __name__ == "__main__":
    main()
