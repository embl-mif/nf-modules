# ZARRSMITH_ENDOMORPH_MERGE

Merges the per-dataset endo_morph CSVs into one per-cell CSV, one per-image summary CSV and one per-gap CSV (header of the first file kept), and copies them to `merged_dir`.

## Input
| Name | Type | Description |
| --- | --- | --- |
| `cells_csvs` | path | All per-cell CSVs |
| `summary_csvs` | path | All per-image summary CSVs |
| `gaps_csvs` | path | All per-gap CSVs |
| `merged_dir` | val | Directory the merged files are published to |

## Output
| Name | Type | Description |
| --- | --- | --- |
| (unnamed) | path | `endomorph_cells_merged.csv` |
| (unnamed) | path | `endomorph_summary_merged.csv` |
| (unnamed) | path | `endomorph_gaps_merged.csv` |
