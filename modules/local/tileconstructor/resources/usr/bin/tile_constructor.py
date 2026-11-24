#!/usr/bin/env python3

"""Cut an OME-Zarr image into regular tiles and write their bounds to a CSV.

The tiles are planned by ZarrSmith (`zarrsmith.utils.tiling.plan_tiles`), from the
image as ZarrSmith reads it (shape, scale and NGFF translation of the chosen
resolution level), so the tiles agree with the later ZarrSmith steps
(ZARRSMITH_FUSE_LABELS).

Tile sizes: an axis given a --tile-size-<ax> is cut at that size. An axis without
one is cut only when a tile would read more than --max-tile-voxels voxels (halo
included, times --resize ** ndim): ZarrSmith then picks the tile counts with the least
cut area that fit. Without --max-tile-voxels, an axis without a tile size is one tile.

Bounds are in voxels of the chosen level, as half-open intervals [min, max). Each row
has a unique dataset_id, <source_dataset_id>__<tile_id>. Tiles are numbered from 1, x
fastest, then y, then z.
"""

import argparse
import csv
import math

from zarrsmith.microscope_manager import ZarrManager
from zarrsmith.utils.tiling import plan_tiles

FIELDNAMES = [
    "dataset_id", "source_dataset_id", "tile_id", "tile_index", "input_uri",
    "x_min", "x_max", "y_min", "y_max", "z_min", "z_max", "resolution_level",
    "tile_overlap_x", "tile_overlap_y", "tile_overlap_z",
    "tile_size_x", "tile_size_y", "tile_size_z",
]


def _level_geometry(image: str, level: int) -> tuple[dict, dict, dict]:
    """(shape, scale, translation) of the spatial axes of ``level``, by axis name."""
    md = ZarrManager(image, mode="r").metadata
    n_levels = len(md["size"])
    if not 0 <= level < n_levels:
        raise ValueError(f"--resolution-level {level} is out of range for an input with {n_levels} level(s).")
    axes = list(md["axes"])
    spatial = [ax for ax in axes if ax in "zyx"]
    if "x" not in spatial or "y" not in spatial:
        raise ValueError(f"The input must have x and y axes, got {''.join(axes)!r}.")
    shape = dict(zip(axes, md["size"][level]))
    scale = dict(zip(spatial, md["scales"][level]))
    translations = md.get("translations") or []
    translation = dict(translations[level]) if level < len(translations) and translations[level] else {}
    return {ax: int(shape[ax]) for ax in spatial}, scale, translation


def _roi(args, axes, physical: bool) -> dict:
    """``{axis: (start, stop)}`` from --<ax>-min/--<ax>-max; voxel bounds widened to whole voxels."""
    roi = {}
    for ax in "zyx":
        lo, hi = getattr(args, f"{ax}_min"), getattr(args, f"{ax}_max")
        if lo is None and hi is None:
            continue
        if ax not in axes:
            raise ValueError(f"The input has no {ax} axis, but --{ax}-min/--{ax}-max was given.")
        if lo is not None and hi is not None and hi <= lo:
            raise ValueError(f"{ax} bounds are invalid: min={lo} must be < max={hi}.")
        if not physical:
            lo = None if lo is None else int(math.floor(lo + 1e-9))
            hi = None if hi is None else int(math.ceil(hi - 1e-9))
        roi[ax] = (lo, hi)
    return roi


def _whole(value, name):
    """Voxel lengths must be whole numbers for plan_tiles; accept 512.0 from Nextflow."""
    if value is None:
        return None
    if float(value) != int(value):
        raise ValueError(f"{name} in voxels must be a whole number, got {value}; use --use-physical-units for units.")
    return int(value)


def construct_tiles(args) -> None:
    physical = bool(args.use_physical_units)
    shape, scale, translation = _level_geometry(args.input_image_path, args.resolution_level)
    axes = tuple(ax for ax in "zyx" if ax in shape)

    sizes = {ax: getattr(args, f"tile_size_{ax}") for ax in axes}
    sizes = {ax: v for ax, v in sizes.items() if v is not None}
    overlap, halo = args.tile_overlap, args.halo
    if not physical:
        sizes = {ax: _whole(v, f"--tile-size-{ax}") for ax, v in sizes.items()}
        overlap, halo = _whole(overlap, "--tile-overlap"), _whole(halo, "--halo")

    plan = plan_tiles(
        shape,
        tile_size=sizes or None,
        overlap=overlap,
        halo=halo,
        roi=_roi(args, axes, physical) or None,
        scale=scale,
        translation=translation,
        physical=physical,
        level=args.resolution_level,
        max_tile_voxels=args.max_tile_voxels,
        resize=args.resize,
    )

    rows = []
    for tile in plan.tiles:
        row = {key: "" for key in FIELDNAMES}
        row.update({
            "dataset_id": f"{args.source_dataset_id}__{tile.tile_id}",
            "source_dataset_id": args.source_dataset_id,
            "tile_id": tile.tile_id,
            "tile_index": tile.index,
            "input_uri": args.input_image_path,
            "resolution_level": plan.level,
        })
        for ax in plan.axes:
            row[f"{ax}_min"], row[f"{ax}_max"] = tile.bounds[ax]
            row[f"tile_overlap_{ax}"] = plan.overlap[ax]
            # the planned tile size (constant over the rows of one image), not this tile's own
            row[f"tile_size_{ax}"] = plan.tile_size[ax]
        rows.append(row)

    with open(args.output_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    counts = " x ".join(f"{plan.counts[ax]} ({ax})" for ax in plan.axes)
    print(f"Wrote {len(rows)} tile(s) to {args.output_csv}: {counts}, tile size {plan.tile_size}, "
          f"overlap {plan.overlap}, halo {plan.halo}, largest read {plan.max_load_voxels} voxels")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--source-dataset-id", required=True)
    parser.add_argument("--input-image-path", required=True)
    parser.add_argument("--output-csv", required=True)

    for ax in "xyz":
        parser.add_argument(f"--tile-size-{ax}", type=float, default=None,
                            help=f"Tile size along {ax}. Default: not split, unless --max-tile-voxels cuts it.")
    parser.add_argument("--max-tile-voxels", type=int, default=None,
                        help="Most voxels one tile may read, halo included (e.g. what fits in GPU or host "
                             "memory). Axes without a tile size are cut to fit it.")
    parser.add_argument("--resize", type=float, default=1.0,
                        help="Linear factor by which the tool enlarges each tile before processing "
                             "(e.g. 30 / diameter for Cellpose); the budget is checked on voxels * resize ** ndim.")
    parser.add_argument("--tile-overlap", type=float, default=0.0)
    parser.add_argument("--resolution-level", type=int, default=0)
    parser.add_argument("--halo", type=float, default=0.0,
                        help="Context read around each tile by the next step (same units as --tile-overlap). "
                             "Counts towards --max-tile-voxels, and keeps a last tile from being smaller than it.")

    for ax in "xyz":
        parser.add_argument(f"--{ax}-min", type=float, default=None)
        parser.add_argument(f"--{ax}-max", type=float, default=None)

    parser.add_argument("--use-physical-units", action="store_true",
                        help="Tile sizes, overlap, halo and ROI bounds are in physical units "
                             "(the NGFF translation is the centre of the first voxel).")
    return parser


def main() -> None:
    construct_tiles(build_arg_parser().parse_args())


if __name__ == "__main__":
    main()
