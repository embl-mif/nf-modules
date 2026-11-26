process TILE_CONSTRUCTOR {
    tag "$meta.id"
    label 'process_low'

    // `zarrsmith plan-tiles`, from the conda environment in environment.yml (a commit of ZarrSmith's
    // branch consolidated_translation, not released yet; pin a release once it is out).
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(local_image), val(image_uri)

    output:
    tuple val(meta), path("*_tiles.csv"), emit: tiles

    script:
    def image = image_uri ?: local_image
    def output_csv = "${meta.id}_tiles.csv"
    def shellQuote = { value -> "'${value.toString().replace("'", "'\"'\"'")}'" }
    def present = { value -> value != null && value.toString().trim() }
    // tile sizes per axis ("x=512;y=512"), and the part of the level to tile ("x=0:900;z=10:")
    def tileSize = ['x', 'y', 'z'].findAll { present(meta["tile_size_${it}"]) }
        .collect { "${it}=${meta["tile_size_${it}"]}" }.join(';')
    def roi = ['x', 'y', 'z'].findAll { present(meta["${it}_min"]) || present(meta["${it}_max"]) }
        .collect { "${it}=${present(meta["${it}_min"]) ? meta["${it}_min"] : ''}:${present(meta["${it}_max"]) ? meta["${it}_max"] : ''}" }
        .join(';')
    def optional = { flag, value -> present(value) ? "${flag} ${shellQuote(value)}" : "" }

    // zarrsmith plan-tiles: one row per tile, with its bounds and its forward neighbours
    // (right/bottom/behind_dataset_id), in the voxels of --level unless --physical.
    """
    zarrsmith plan-tiles \\
        -i ${shellQuote(image)} \\
        -oc ${output_csv} \\
        -id ${shellQuote(meta.id)} \\
        ${optional('-lv', meta.resolution_level)} \\
        ${optional('-ts', tileSize)} \\
        ${optional('-mtv', meta.max_tile_voxels)} \\
        ${optional('-rs', meta.tile_resize)} \\
        ${optional('-ov', meta.tile_overlap)} \\
        ${optional('-ha', meta.halo)} \\
        ${optional('-roi', roi)} \\
        ${meta.use_physical_units ? '-pu' : ''}
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
