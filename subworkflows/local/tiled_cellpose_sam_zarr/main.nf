include { TILE_CONSTRUCTOR      } from '../../../modules/local/tileconstructor/main'
include { CELLPOSE_SAM_ZARR_MIF } from '../../../modules/local/cellpose_sam_zarr_mif/main'
include { LOCAL_LABEL_MATCHER   } from '../../../modules/local/local_label_matcher/main'
include { LABEL_SOLVER          } from '../../../modules/local/label_solver/main'
include { ZARRSMITH_FUSE_LABELS } from '../../../modules/local/zarrsmith_fuse_labels/main'

workflow TILED_CELLPOSE_SAM_ZARR {

    take:
    ch_images
    val_label_offset_step
    val_multicut_beta

    main:

    ch_versions = Channel.empty()

    def repoRootDir = workflow.projectDir.getParent().getParent()
    def noNeighborZarr  = file("${repoRootDir}/assets/NO_FILE")
    def noImageSentinel = file("${repoRootDir}/assets/NO_FILE")
    def tileConstructorKeys = ['halo', 'resolution_level', 'tile_overlap', 'tile_size_x', 'tile_size_y', 'tile_size_z',
                               'use_physical_units', 'x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max']
    def cellposeKeys = ['diameter', 'niter', 'segmentation_channels', 'flow_threshold', 'cellprob_threshold', 'halo']

    ch_split_images = ch_images.map { meta, image ->
        def (localImage, imageUri) = Utils.splitImageInput(image, noImageSentinel)
        tuple(meta.id, meta, image, localImage, imageUri)
    }

    // Tile budget: the most voxels one Cellpose tile may read, halo included. Axes without a
    // tile_size_<ax> are cut by ZarrSmith's planner to fit it. The default fits the first
    // attempt of the 'gpu' label in conf/embl_hd.config (about 2.7 GB + 329 bytes per input
    // voxel, measured at diameter 15, so per input voxel: tile_resize stays 1 by default).
    // Set it per image (CSV column max_tile_voxels) or for the run (--max_tile_voxels).
    // tile_resize (e.g. 30 / diameter) multiplies the voxels counted against the budget.
    def defaultMaxTileVoxels = 300_000_000L
    def tileConstructorMeta = { meta ->
        def budget = Utils.normalizeOptionalText(meta.max_tile_voxels) ?: Utils.normalizeOptionalText(params.get('max_tile_voxels'))
        def resize = Utils.normalizeOptionalText(meta.tile_resize)
        meta.subMap(tileConstructorKeys) + [
            id             : meta.id,
            max_tile_voxels: budget ? new BigDecimal(budget).toBigInteger() : defaultMaxTileVoxels,
            tile_resize    : resize ? resize as double : null,
        ]
    }

    TILE_CONSTRUCTOR( ch_split_images.map { id, meta, image, localImage, imageUri -> tuple(tileConstructorMeta(meta), localImage, imageUri) } )

    ch_image_tiles = TILE_CONSTRUCTOR.out.tiles
        .map { parentMeta, tileCsv -> tuple(parentMeta.id, tileCsv) }
        .join( ch_split_images.map { id, meta, image, localImage, imageUri -> tuple(id, meta, image, localImage, imageUri) } )
        .flatMap { _srcId, tileCsv, parentMeta, image, localImage, imageUri ->
            tileCsv.splitCsv(header: true).collect { row ->
                def tileMeta = parentMeta.subMap(cellposeKeys) + [
                    id: row.dataset_id,
                    source_dataset_id: row.source_dataset_id,
                    tile_id: row.tile_id,
                    tile_index: row.tile_index,
                    resolution_level: row.resolution_level,
                    x_min: row.x_min,
                    x_max: row.x_max,
                    y_min: row.y_min,
                    y_max: row.y_max,
                    z_min: row.z_min,
                    z_max: row.z_max,
                    tile_overlap_x: row.tile_overlap_x,
                    tile_overlap_y: row.tile_overlap_y,
                    tile_overlap_z: row.tile_overlap_z,
                    tile_size_x: row.tile_size_x,
                    tile_size_y: row.tile_size_y,
                    tile_size_z: row.tile_size_z,
                    // forward neighbours from the tile plan (blank at the end of an axis)
                    right_dataset_id : row.right_dataset_id,
                    right_tile_index : row.right_tile_index,
                    bottom_dataset_id: row.bottom_dataset_id,
                    bottom_tile_index: row.bottom_tile_index,
                    behind_dataset_id: row.behind_dataset_id,
                    behind_tile_index: row.behind_tile_index,
                    use_physical_units: false,
                    do_3D: parentMeta.do_3D ?: false,
                    zarr_output: "${row.dataset_id}_labels.ome.zarr",
                    is_tiled: true,
                ]

                tuple(tileMeta, image, localImage, imageUri)
            }
        }

    ch_cellpose_model = params.cellpose_model_path
        ? Channel.value(file(params.cellpose_model_path, checkIfExists: true))
        : Channel.value(file("${repoRootDir}/assets/NO_MODEL_FILE"))

    CELLPOSE_SAM_ZARR_MIF(
        ch_image_tiles.map { meta, image, localImage, imageUri -> tuple(meta, localImage, imageUri) },
        ch_cellpose_model
    )

    ch_all_labels = CELLPOSE_SAM_ZARR_MIF.out.labels

    ch_tiles_and_segmentations_csv = ch_all_labels
        .map { meta, label ->
            tuple(meta.source_dataset_id, "${meta.id},${meta.source_dataset_id},${meta.tile_id},${meta.tile_index},${meta.x_min},${meta.x_max},${meta.y_min},${meta.y_max},${meta.z_min},${meta.z_max},${meta.resolution_level},${meta.tile_overlap_x},${meta.tile_overlap_y},${meta.tile_overlap_z},${label}")
        }
        .collectFile(seed: 'dataset_id,source_dataset_id,tile_id,tile_index,x_min,x_max,y_min,y_max,z_min,z_max,resolution_level,tile_overlap_x,tile_overlap_y,tile_overlap_z,segmented_tile_path', newLine: true) { srcId, row ->
            [ "${srcId}_tiles_and_segmentations.csv", row ]
        }
        .map { csvPath -> tuple(csvPath.baseName.replace('_tiles_and_segmentations', ''), csvPath) }

    ch_label_path_by_id = ch_all_labels
        .map { meta, label -> tuple(meta.id, label) }
        .toList()
        .map { pairs -> pairs.collectEntries { id, label -> [(id): label] } }

    // Every tile's meta already names its forward neighbours (TILE_CONSTRUCTOR, from the
    // tile plan); here they are only resolved to the neighbours' segmented stores.
    ch_local_matcher_input = ch_all_labels
        .combine(ch_label_path_by_id)
        .map { meta, selfLabel, labelPathById ->
            def rightZarr  = labelPathById[meta.right_dataset_id]  ?: noNeighborZarr
            def bottomZarr = labelPathById[meta.bottom_dataset_id] ?: noNeighborZarr
            def behindZarr = labelPathById[meta.behind_dataset_id] ?: noNeighborZarr

            tuple(meta + [label_offset_step: val_label_offset_step], selfLabel, rightZarr, bottomZarr, behindZarr)
        }

    LOCAL_LABEL_MATCHER(ch_local_matcher_input)


    ch_local_graphs_grouped = LOCAL_LABEL_MATCHER.out.local_graph
        .map { meta, json -> tuple(meta.source_dataset_id, json) }
        .groupTuple()
        .map { srcId, jsonList -> tuple([id: srcId, beta: val_multicut_beta], jsonList) }

    LABEL_SOLVER(ch_local_graphs_grouped)

    // One fusion per source image: every tile store (placed by its own NGFF translation),
    // the tile indices that set each tile's label offset, and the solver's mapping.
    ch_fuse_input = ch_all_labels
        .map { meta, label -> tuple(meta.source_dataset_id, label) }
        .groupTuple()
        .join( ch_tiles_and_segmentations_csv )
        .join( LABEL_SOLVER.out.global_mapping.map { meta, json -> tuple(meta.id, json) } )
        .join( ch_split_images.map { id, meta, image, localImage, imageUri -> tuple(id, meta, localImage, imageUri) } )
        .map { srcId, labels, tilesCsv, mapping, parentMeta, localImage, imageUri ->
            def mergedMeta = [
                id               : srcId,
                dtype            : 'uint32',
                data_type        : 'label',
                zarr_output      : Utils.resolveOutput(parentMeta, 'labels_output', "${srcId}_merged_labels.ome.zarr", srcId),
                label_offset_step: val_label_offset_step,
                resolution_level : parentMeta.resolution_level,
            ]
            tuple(mergedMeta, labels, tilesCsv, mapping, localImage, imageUri)
        }

    ZARRSMITH_FUSE_LABELS(ch_fuse_input)
    ch_versions = ch_versions.mix(ZARRSMITH_FUSE_LABELS.out.versions)

    emit:
    labels   = ZARRSMITH_FUSE_LABELS.out.zarr
    versions = ch_versions
}
