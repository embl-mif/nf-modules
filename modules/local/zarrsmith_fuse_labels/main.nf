process ZARRSMITH_FUSE_LABELS {
    tag "${meta.id}"
    label 'process_medium'
    // `zarrsmith fuse-labels`, from the conda environment in environment.yml (a commit of ZarrSmith's
    // branch consolidated_translation, not released yet; pin a release once it is out).
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(tile_labels, stageAs: 'tiles/*'), path(tiles_csv), path(global_mapping), path(local_image), val(image_uri)

    output:
    tuple val(meta), path(output_zarr_link), emit: zarr
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    output_zarr = meta.zarr_output ?: "${meta.id}_merged_labels.ome.zarr"
    // A `path` output must live in the work dir: an absolute output is exposed there as a
    // symlink named after it (see Utils.shellLinkOutput), a relative one is the output itself.
    output_zarr_link = new File(output_zarr.toString()).name

    if (meta.label_offset_step == null) {
        error "ZARRSMITH_FUSE_LABELS (${meta.id}): meta.label_offset_step is required."
    }
    def args = task.ext.args ?: ''
    def shellQuote = { value -> "'${value.toString().replace("'", "'\"'\"'")}'" }
    // local_image is a NO_*-named sentinel when the image is given by URI (see PREALLOCATE_ZARR)
    def source = image_uri ?: (!local_image.name.startsWith('NO_') ? local_image : null)

    // -i: a folder of OME-Zarr stores, each placed by its own NGFF translation.
    // --tiles_csv/--label_offset_step: every tile's ids get tile_index * step added (the numbering
    // LOCAL_LABEL_MATCHER gave the graph nodes), after checking that the staged stores and the CSV
    // list the same tiles; the offsets used are written next to the output (<name>.label_offsets.csv).
    // -lm: LABEL_SOLVER's mapping merges the objects split between tiles. Labels are never blended.
    // -idf/-sl: the pyramid follows the source image from the segmented level on, so level k of the
    // labels matches level resolution_level + k of the image; -nl/-df override it.
    def fuseParts = [
        'zarrsmith', 'fuse-labels',
        '-i', 'tiles',
        '-z', shellQuote(output_zarr),
        '--tiles_csv', shellQuote(tiles_csv),
        '--label_offset_step', meta.label_offset_step,
        '-lm', shellQuote(global_mapping),
        '-dt', meta.dtype ?: 'uint32',
        '-rw', task.cpus ?: 1,
        '-ow',
    ]
    if (source)                        { fuseParts.addAll(['-idf', shellQuote(source)]) }
    if (source && meta.resolution_level != null) { fuseParts.addAll(['-sl', meta.resolution_level]) }
    if (meta.num_levels != null)       { fuseParts.addAll(['-nl', meta.num_levels]) }
    if (meta.downscale_factor != null) { fuseParts.add('-df'); fuseParts.addAll(meta.downscale_factor.toString().split(/[\s,;]+/).findAll { it }) }
    if (meta.pyramid_method != null)   { fuseParts.addAll(['-me', meta.pyramid_method]) }

    // Keep this script indented with spaces only: Nextflow dedents it, and stray tabs
    // leave the END_VERSIONS heredoc terminator indented, so bash never closes it.
    """
    mkdir -p "\$(dirname ${shellQuote(output_zarr)})"
    ${fuseParts.join(' ')} ${args}

    if [ ! -e ${shellQuote(output_zarr)}/zarr.json ] && [ ! -e ${shellQuote(output_zarr)}/.zgroup ]; then
        echo "ERROR: expected zarr store not found: ${output_zarr}" >&2
        exit 1
    fi

    ${Utils.shellLinkOutput(output_zarr)}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        zarrsmith: \$(python -c 'from importlib.metadata import version; print(version("zarrsmith"))' 2>/dev/null || echo "unknown")
    END_VERSIONS
    """

    stub:
    output_zarr = meta.zarr_output ?: "${meta.id}_merged_labels.ome.zarr"
    output_zarr_link = new File(output_zarr.toString()).name
    """
    mkdir -p "${output_zarr}"
    ${Utils.shellLinkOutput(output_zarr)}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        zarrsmith: "stub"
    END_VERSIONS
    """
}
