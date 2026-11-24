#!/usr/bin/env python3

"""Prepare the inputs of `zarrsmith fuse` for ZARRSMITH_FUSE_LABELS.

The fusion itself is done by `zarrsmith fuse --label_offsets ... --label_mapping ...`:
every tile store is placed at its own NGFF translation (written by
CELLPOSE_SAM_ZARR_MIF), its non-zero ids get its offset added, then LABEL_SOLVER's
global mapping merges the objects split between tiles, without blending.

This script writes what the CLI cannot work out on its own:

* the offsets CSV (`store,offset`): offset = tile_index * label_offset_step, the
  numbering LOCAL_LABEL_MATCHER gave the graph nodes. It also checks that the staged
  tile stores and the tiles CSV list the same tiles.
* the pyramid options (`-nl N -df z y x`), written to --pyramid-args: they follow the
  source image from the segmented resolution level on (same number of levels below
  it, same factors), so level k of the labels matches level resolution_level + k of
  the image. Without a source image and without --num-levels/--downscale-factor, the
  file is empty and zarrsmith's defaults apply.
"""

import argparse
import csv
import sys
from pathlib import Path


def _source_pyramid(source: str, level: int) -> tuple[int, tuple[int, int, int]]:
    """Number of levels from ``level`` on, and the (z, y, x) factor between levels."""
    from zarrsmith.microscope_manager import ZarrManager

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
    parser.add_argument("--label-offset-step", type=int, required=True)
    parser.add_argument("--offsets-csv", required=True, help="Output: store,offset CSV for --label_offsets.")
    parser.add_argument("--pyramid-args", required=True, help="Output: -nl/-df options for zarrsmith fuse.")
    parser.add_argument("--source-image", default=None, help="Source image, for the pyramid levels.")
    parser.add_argument("--resolution-level", type=int, default=0, help="Source level that was segmented.")
    parser.add_argument("--num-levels", type=int, default=None, help="Default: as the source image.")
    parser.add_argument("--downscale-factor", type=int, nargs="+", default=None,
                        help="One value or z y x. Default: as the source image.")
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    step = args.label_offset_step
    if step < 1:
        raise ValueError(f"--label-offset-step must be >= 1, got {step}.")

    with open(args.tiles_csv, newline="", encoding="utf-8") as handle:
        index_of = {Path(row["segmented_tile_path"]).name: int(row["tile_index"]) for row in csv.DictReader(handle)}

    stores = sorted(p.name for p in Path(args.tiles_dir).iterdir() if p.is_dir())
    missing = sorted(set(index_of) - set(stores))
    unknown = sorted(set(stores) - set(index_of))
    if missing or unknown:
        raise ValueError(f"Tile stores and {args.tiles_csv} disagree: missing stores {missing}, "
                         f"stores not in the CSV {unknown}.")

    with open(args.offsets_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["store", "offset"])
        for store in stores:
            writer.writerow([store, index_of[store] * step])
    print(f"Wrote offsets of {len(stores)} tile store(s) to {args.offsets_csv} (step {step})", file=sys.stderr)

    levels, factor = (_source_pyramid(args.source_image, args.resolution_level)
                      if args.source_image else (None, None))
    if args.num_levels is not None:
        levels = args.num_levels
    if args.downscale_factor is not None:
        factor = tuple(args.downscale_factor)

    options = []
    if levels is not None:
        options += ["-nl", str(levels)]
    if factor is not None:
        options += ["-df", *(str(f) for f in factor)]
    with open(args.pyramid_args, "w", encoding="utf-8") as handle:
        handle.write(" ".join(options) + "\n")
    print(f"Pyramid options: {' '.join(options) or '(zarrsmith defaults)'}", file=sys.stderr)


if __name__ == "__main__":
    main()
