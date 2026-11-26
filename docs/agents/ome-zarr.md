# OME-Zarr preferences

Written by `/setup-ome-zarr-skills` on 2026-11-26. Edit by hand freely.
Skills read this file when present; when a field is `unknown` they infer from context.

## Spec version

- target: 0.5
- reason: what the pipelines write today. ZarrSmith modules write 0.5 / zarr v3 by default; CELLPOSE_SAM_ZARR_MIF writes its labels in the version of its input; nothing writes 0.6 yet.
- also_read: 0.4 (ilastik modules read and write only zarr v2, so their stores are 0.4 / v2 without shards; `zarr_format 2` is still offered by the conversion workflows), 0.6 (read through ZarrSmith)
- spec: https://ngff.openmicroscopy.org/specifications/0.5/index.html

## RFCs in use

None beyond the target version. RFC-2 (zarr v3) comes with 0.5. RFC-5 (coordinate systems, part of 0.6) is not followed as such, but translations follow its pixel-centre convention through ZarrSmith (see Notes).

## Domain

- modality: light-sheet, confocal/widefield and HCS (Opera) acquisitions converted with ZarrSmith; organoid segmentation (Cellpose-SAM, tiled), ilastik pixel classification, Ultrack cell tracking
- axes: any subset of `tczyx`, in that order
- units: micrometer (space), second (time) (inferred, from ZarrSmith)
- dtypes: uint16 raw (inferred, ZarrSmith default); uint32 labels (merged Cellpose labels); float32 ilastik probabilities (inferred)
- typical_size: unknown (varies by user and dataset); one Cellpose tile reads at most ~300 M voxels incl. halo (max_tile_voxels default)
- extras: omero (from ZarrSmith); labels written as standalone label stores (`*_labels.ome.zarr`, `*_merged_labels.ome.zarr`), not a `labels/` group (inferred); no plate/well, no bioformats2raw.layout
- storage: local cluster filesystems (`/scratch`, `/g/...`); fused labels use ZarrSmith defaults (~4 MB chunks, 256 MB shards on v3); tile planner inner chunk 16 x 256 x 256

## Notes

- Never mix pairs: zarr v2 with NGFF 0.4, zarr v3 with NGFF 0.5.
- Tiles are planned and fused with ZarrSmith (`zarrsmith plan-tiles`, `zarrsmith fuse-labels`); each tile store is placed by its own NGFF translation, read as the centre of the first voxel (0.6 convention).
- Label pyramids are sampled (nearest), never averaged; the merged label pyramid follows the source image from the segmented resolution level on.
- Follow ZarrSmith's own conventions (zarrsmith/docs/agents/ome-zarr.md) for layout, compression and metadata writing.
