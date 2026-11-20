#!/usr/bin/env python3

"""Merge the label tiles of one tiled segmentation into one OME-Zarr label image
(ZARRSMITH_FUSE_LABELS).

Every tile store is read at its NGFF translation (written by CELLPOSE_SAM_ZARR_MIF),
so its place in the merged image comes from its own metadata. Its labels are made
unique the same way LOCAL_LABEL_MATCHER numbered the graph nodes (every non-zero id
+ tile_index * label_offset_step), then LABEL_SOLVER's global mapping merges the
objects split between tiles. ZarrSmith fuses the tiles without blending
(blend="none"): every voxel is copied from one tile, overlaps are cut in their
middle, and the pyramid levels sample voxels instead of averaging them.

The pyramid follows the source image from the segmented resolution level on (same
number of levels below it and the same factors), so level k of the labels matches
level resolution_level + k of the image.
"""

import argparse
import csv
import json
from pathlib import Path

from zarrsmith.microscope_manager import ZarrManager
from zarrsmith.mosaic_manager import RelabeledArray, fuse_tiles, store_tiles


def _source_pyramid(source: str, level: int) -> tuple[int, tuple[int, int, int]]:
    """Number of levels from ``level`` on, and the (z, y, x) factor between levels."""
    md = ZarrManager(source, mode="r").metadata
    scales = md["scales"]
    spatial = [ax for ax in md["axes"] if ax in "zyx"]
    levels = len(scales) - level
    if levels < 1:
        raise ValueError(f"--resolution-level {level} is out of range: the source has {len(scales)} level(s).")
    if levels == 1:
        return 1, (1, 2, 2)
    ratio = {ax: scales[level + 1][i] / scales[level][i] for i, ax in enumerate(spatial)}
    return levels, tuple(max(1, round(ratio.get(ax, 1.0))) for ax in "zyx")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--tiles-dir", required=True, help="Folder holding the tile label stores.")
    parser.add_argument("--tiles-csv", required=True,
                        help="tiles_and_segmentations.csv: tile_index and segmented_tile_path of every tile.")
    parser.add_argument("--global-mapping", required=True, help="LABEL_SOLVER's global_label_mapping.json.")
    parser.add_argument("--label-offset-step", type=int, required=True)
    parser.add_argument("--output-zarr", required=True)
    parser.add_argument("--source-image", default=None, help="Source image, for the pyramid levels.")
    parser.add_argument("--resolution-level", type=int, default=0, help="Source level that was segmented.")
    parser.add_argument("--num-levels", type=int, default=None, help="Default: as the source image.")
    parser.add_argument("--downscale-factor", type=int, nargs="+", default=None,
                        help="One value or z y x. Default: as the source image, else 2.")
    parser.add_argument("--dtype", default="uint32")
    parser.add_argument("--name", default=None)
    parser.add_argument("--workers", type=int, default=1)
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    step = args.label_offset_step

    with open(args.tiles_csv, newline="", encoding="utf-8") as handle:
        index_of = {Path(row["segmented_tile_path"]).name: int(row["tile_index"]) for row in csv.DictReader(handle)}
    with open(args.global_mapping, encoding="utf-8") as handle:
        mapping = {int(k): int(v) for k, v in json.load(handle).items()}

    stores = sorted(p for p in Path(args.tiles_dir).iterdir() if p.is_dir())
    missing = sorted(set(index_of) - {p.name for p in stores})
    unknown = sorted({p.name for p in stores} - set(index_of))
    if missing or unknown:
        raise ValueError(f"Tile stores and {args.tiles_csv} disagree: missing stores {missing}, "
                         f"stores not in the CSV {unknown}.")

    tiles = []
    for store in stores:
        store_index = index_of[store.name]
        offset = store_index * step
        # only this tile's node ids, so each worker gets a small mapping
        own = {k: v for k, v in mapping.items() if offset < k < offset + step}
        found, _placed = store_tiles(store)
        for tile in found:
            tile.data = RelabeledArray(tile.data, offset=offset, mapping=own, dtype=args.dtype)
            tiles.append(tile)
        print(f"{store.name}: tile {store_index}, offset {offset}, {len(own)} merged ids")

    levels, factor = (_source_pyramid(args.source_image, args.resolution_level)
                      if args.source_image else (1, (1, 2, 2)))
    if args.num_levels is not None:
        levels = args.num_levels
    if args.downscale_factor is not None:
        factor = args.downscale_factor[0] if len(args.downscale_factor) == 1 else tuple(args.downscale_factor)
    print(f"Fusing {len(tiles)} tile(s) into {args.output_zarr}: {levels} level(s), factor {factor}")

    fuse_tiles(
        tiles, args.output_zarr,
        blend="none", data_type="label", out_dtype=args.dtype,
        levels=levels, downscale_factor=factor,
        workers=args.workers, name=args.name, overwrite=True,
        transforms_path=False, progress=False,
    )
    print("Done.")


if __name__ == "__main__":
    main()
