process CUSTOM_DUMPSOFTWAREVERSIONS {
    label 'process_single'

    // Requires `pyyaml` which does not have a dedicated container but is in the MultiQC container
    conda (params.enable_conda ? 'bioconda::multiqc=1.23' : null)
    container "wslh-bioinformatics/multiqc@sha256:b0380bd93470fa15d30ec537ab9385b504addd0f57a448601215f499e2d43bd1"

    input:
    path versions

    output:
    path "software_versions.yml"    , emit: yml
    path "software_versions_mqc.yml", emit: mqc_yml
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    template 'dumpsoftwareversions.py'
}
