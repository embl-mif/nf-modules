process ZARRSMITH_FUSE_LABELS {
    tag "${meta.id}"
    label 'process_medium'
    // Needs ZarrSmith with `zarrsmith fuse --label_offsets/--label_mapping` (branch
    // consolidated_translation, not released yet). Pin this to the release version once it is out.
    container 'ghcr.io/embl-mif/zarrsmith:consolidated_translation'
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

    // offsets CSV (tile_index * label_offset_step per store) and the source's pyramid (-nl/-df)
    def prepareParts = [
        'fuse_labels_inputs.py',
        '--tiles-dir', 'tiles',
        '--tiles-csv', shellQuote(tiles_csv),
        '--label-offset-step', meta.label_offset_step,
        '--offsets-csv', 'label_offsets.csv',
        '--pyramid-args', 'pyramid_args.txt',
    ]
    if (source)                        { prepareParts.addAll(['--source-image', shellQuote(source)]) }
    if (meta.resolution_level != null) { prepareParts.addAll(['--resolution-level', meta.resolution_level]) }
    if (meta.num_levels != null)       { prepareParts.addAll(['--num-levels', meta.num_levels]) }
    if (meta.downscale_factor != null) { prepareParts.add('--downscale-factor'); prepareParts.addAll(meta.downscale_factor.toString().split(/[\s,;]+/).findAll { it }) }

    // -m zarr: a folder of OME-Zarr stores, each placed by its own NGFF translation.
    // --label_offsets/--label_mapping imply --data_type label and blend none.
    def fuseParts = [
        'zarrsmith', 'fuse',
        '-m', 'zarr',
        '-i', 'tiles',
        '-z', shellQuote(output_zarr),
        '-lo', 'label_offsets.csv',
        '-lm', shellQuote(global_mapping),
        '-od', meta.dtype ?: 'uint32',
        '-rw', task.cpus ?: 1,
        '-ow',
    ]
    if (meta.pyramid_method != null)   { fuseParts.addAll(['-me', meta.pyramid_method]) }

    // Keep this script indented with spaces only: Nextflow dedents it, and stray tabs
    // leave the END_VERSIONS heredoc terminator indented, so bash never closes it.
    """
    mkdir -p "\$(dirname ${shellQuote(output_zarr)})"
    ${prepareParts.join(' ')}
    ${fuseParts.join(' ')} \$(cat pyramid_args.txt) ${args}

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
