process TILE_CONSTRUCTOR {
    tag "$meta.id"
    label 'process_low'

    // Needs ZarrSmith's tile planner (zarrsmith.utils.tiling.plan_tiles, branch
    // consolidated_translation, not released yet). Pin this to the release version once it is out.
    container 'ghcr.io/embl-mif/zarrsmith:consolidated_translation'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(local_image), val(image_uri)

    output:
    tuple val(meta), path("*_tiles.csv"), emit: tiles

    script:
    def image = image_uri ?: local_image
    def output_csv = "${meta.id}_tiles.csv"
    def optional = { flag, value -> value != null ? "${flag} ${value}" : "" }

    """
    tile_constructor.py \\
        --source-dataset-id "${meta.id}" \\
        --input-image-path \"${image}\" \\
        --output-csv ${output_csv} \\
        ${optional('--tile-size-x', meta.tile_size_x)} \\
        ${optional('--tile-size-y', meta.tile_size_y)} \\
        ${optional('--tile-size-z', meta.tile_size_z)} \\
        ${optional('--max-tile-voxels', meta.max_tile_voxels)} \\
        ${optional('--resize', meta.tile_resize)} \\
        ${optional('--tile-overlap', meta.tile_overlap)} \\
        ${optional('--halo', meta.halo)} \\
        ${optional('--resolution-level', meta.resolution_level)} \\
        ${optional('--x-min', meta.x_min)} \\
        ${optional('--x-max', meta.x_max)} \\
        ${optional('--y-min', meta.y_min)} \\
        ${optional('--y-max', meta.y_max)} \\
        ${optional('--z-min', meta.z_min)} \\
        ${optional('--z-max', meta.z_max)} \\
        ${meta.use_physical_units ? '--use-physical-units' : ''}
    """

    stub:
    def image = image_uri ?: local_image
    def output_csv = "${meta.id}_tiles.csv"
    """
    cat > ${output_csv} <<-END_TILE_CSV
    dataset_id,source_dataset_id,tile_id,tile_index,input_uri,x_min,x_max,y_min,y_max,z_min,z_max,resolution_level,tile_overlap_x,tile_overlap_y,tile_overlap_z,tile_size_x,tile_size_y,tile_size_z
    ${meta.id}__tile_000001,${meta.id},tile_000001,1,${image},0,1,0,1,,,0,0,0,,1,1,
    END_TILE_CSV
    """
}
