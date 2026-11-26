process ZARRSMITH_FEATURE_EXTRACTION {
    tag "${meta.id}"
    label 'process_medium'
    // `zarrsmith measure`, from the conda environment in environment.yml (a commit of ZarrSmith's
    // branch consolidated_translation, not released yet; pin a release once it is out).
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(intensity_zarr), path(labels_zarr)

    output:
    tuple val(meta), path(output_csv_link), emit: csv
    path "versions.yml"                   , emit: versions

    script:
    output_csv = meta.csv_output ?: "${meta.id}_features.csv"
    // A `path` output must live in the work dir: an absolute output is exposed there as a
    // symlink named after it (see Utils.shellLinkOutput), a relative one is the output itself.
    output_csv_link = new File(output_csv.toString()).name

    def shellQuote = { value ->
        "'${value.toString().replace("'", "'\"'\"'")}'"
    }

    // List values ("0,2", "solidity; euler_number") become one argument per item.
    def listArg = { flag, value ->
        def items = value ? value.toString().trim().split(/[;,\s]+/).findAll { it } : []
        items ? "${flag} ${items.collect { shellQuote(it) }.join(' ')}" : ""
    }
    def intensityLevelArg = meta.features_intensity_resolution_level != null ? "-lv ${meta.features_intensity_resolution_level}" : ""
    def labelsLevelArg = meta.features_labels_resolution_level != null ? "-llv ${meta.features_labels_resolution_level}" : ""
    def channelsArg = listArg('-ch', meta.features_channels)
    def extraPropertiesArg = listArg('--properties', meta.features_extra_properties)
    def extraIntensityPropertiesArg = listArg('--intensity_properties', meta.features_extra_intensity_properties)

    """
    if [ -e ${shellQuote(output_csv)} ]; then
        echo "Output already exists, skipping feature extraction: ${output_csv}"
    else
        mkdir -p \$(dirname ${shellQuote(output_csv)})

        zarrsmith measure \
            -i "${intensity_zarr}" \
            -l "${labels_zarr}" \
            -oc "${output_csv}" \
            --dataset_id ${shellQuote(meta.id)} \
            ${intensityLevelArg} \
            ${labelsLevelArg} \
            ${channelsArg} \
            ${extraPropertiesArg} \
            ${extraIntensityPropertiesArg}
    fi

    if [ ! -e ${shellQuote(output_csv)} ]; then
        echo "ERROR: expected features CSV not found: ${output_csv}" >&2
        exit 1
    fi

    ${Utils.shellLinkOutput(output_csv)}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        zarrsmith: \$(python -c 'from importlib.metadata import version; print(version("zarrsmith"))' 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
