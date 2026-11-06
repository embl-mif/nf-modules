# ZARRSMITH_RESAMPLE

Changes the pixel size of an OME-Zarr image with `zarrsmith resample`, rescaling the data accordingly. The copy keeps the compression and number of levels of the input, with chunks and shards fitted to each level. Image subgroups and labels are resampled too, unless turned off.

Images use linear interpolation, with anti-aliasing on the axes that shrink; labels use nearest neighbour. The pixel size is in the units of the dataset, usually micrometer.

Requires ZarrSmith 0.6.5 or later.

If `zarr_output` already holds a complete OME-Zarr, the step is skipped; an incomplete one is removed and written again.

## Input
| Name | Type | Description |
| --- | --- | --- |
| `meta` | val | Sample map, see below |
| `input_zarr` | path | Input OME-Zarr |

### `meta` keys
| Key | Required | ZarrSmith option | Description |
| --- | --- | --- | --- |
| `id` | yes | | Sample name |
| `pixel_size` | yes | `--pixel_size` | New pixel size: one value for every spatial axis (`1`), or three in Z Y X order (`1 0.5 0.5`; spaces, commas or semicolons) |
| `zarr_output` | no | `--zarr_path` | Output path. Default: `<id>_resampled.ome.zarr` in the work directory |
| `fast` | no | `--fast` | `true` to read images from the coarsest pyramid level still fine enough instead of level 0: much less data when shrinking, but keeps the aliasing of that level. Default: `false` |
| `include_groups` | no | `--include_groups` | `false` to skip image subgroups. Default: `true` |
| `include_labels` | no | `--include_labels` | `false` to skip labels. Default: `true` |
| `num_workers` | no | `--num_workers` | Writer threads. Default: ZarrSmith picks (half the cores, at most 8, fewer if memory requires it) |

Extra options can be passed with `ext.args`.

## Output
| Name | Type | Description |
| --- | --- | --- |
| `zarr` | `tuple(meta, zarr)` | Resampled OME-Zarr (`meta.zarr_output`) |
| `versions` | path | `versions.yml` |

## Example
```nextflow
include { ZARRSMITH_RESAMPLE } from '../../modules/local/zarrsmith_resample/main'

ch = Channel.of(
    [[id: 'brain', pixel_size: '1', zarr_output: '/scratch/brain_1um.ome.zarr'], file('/data/brain.ome.zarr')]
)
ZARRSMITH_RESAMPLE(ch)
```
