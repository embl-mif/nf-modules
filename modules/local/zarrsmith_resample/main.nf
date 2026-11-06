process ZARRSMITH_RESAMPLE {
    tag "${meta.id}"
    label 'process_medium'
    container 'ghcr.io/embl-mif/zarrsmith:0.6.5'

    input:
    tuple val(meta), path(input_zarr)

    output:
    tuple val(meta), path(output_zarr_link), emit: zarr
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    output_zarr = meta.zarr_output ?: "${meta.id}_resampled.ome.zarr"
    // A `path` output must live in the work dir: an absolute output is exposed there as a
    // symlink named after it (see Utils.shellLinkOutput), a relative one is the output itself.
    output_zarr_link = new File(output_zarr.toString()).name

    def args = task.ext.args ?: ''

    def isPresent = { value ->
        value != null &&
        value.toString().trim() &&
        !['-1', 'none', 'null', 'nan'].contains(value.toString().trim().toLowerCase())
    }
    def isTrue = { value ->
        isPresent(value) && ['true', '1', 'yes', 'y'].contains(value.toString().trim().toLowerCase())
    }
    def shellQuote = { value ->
        "'${value.toString().replace("'", "'\"'\"'")}'"
    }

    // pixel_size: one value for every spatial axis, or three in Z Y X order,
    // separated by spaces, commas or semicolons, in the units of the dataset.
    if (!isPresent(meta.pixel_size)) {
        error "ZARRSMITH_RESAMPLE (${meta.id}): meta.pixel_size is required (one value, or three in Z Y X order)."
    }
    def pixelSizes = meta.pixel_size.toString().trim().split(/[\s,;]+/).findAll { it }
    if (!(pixelSizes.size() in [1, 3])) {
        error "ZARRSMITH_RESAMPLE (${meta.id}): invalid pixel_size '${meta.pixel_size}': give one value, or three in Z Y X order."
    }

    def commandParts = [
        'zarrsmith', 'resample',
        '--input_path', shellQuote(input_zarr),
        '--zarr_path', shellQuote(output_zarr),
        '--pixel_size'
    ] + pixelSizes.collect { shellQuote(it) }

    if (isTrue(meta.fast)) {
        commandParts.add('--fast')
    }
    [
        include_groups: '--include_groups',
        include_labels: '--include_labels'
    ].each { key, option ->
        if (isPresent(meta[key])) {
            commandParts.addAll([option, shellQuote(meta[key].toString().trim().toLowerCase())])
        }
    }
    if (isPresent(meta.num_workers)) {
        commandParts.addAll(['--num_workers', shellQuote(meta.num_workers)])
    }

    def command = commandParts.join(' ')

    // Keep this script indented with spaces only: Nextflow dedents it, and stray tabs
    // leave the END_VERSIONS heredoc terminator indented, so bash never closes it.
    """
    out=${shellQuote(output_zarr)}
    if [ -d "\$out" ] && { grep -q '"ome"' "\$out/zarr.json" 2>/dev/null || grep -q '"multiscales"' "\$out/.zattrs" 2>/dev/null; }; then
        echo "Output already exists, skipping resampling: ${output_zarr}"
    else
        if [ -e "\$out" ]; then
            echo "Incomplete output found (no OME metadata), removing: ${output_zarr}"
            rm -rf "\$out"
        fi
        mkdir -p "\$(dirname "\$out")"
        ${command} ${args}
    fi

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
    output_zarr = meta.zarr_output ?: "${meta.id}_resampled.ome.zarr"
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
