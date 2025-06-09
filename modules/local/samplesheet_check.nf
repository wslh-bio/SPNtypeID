process SAMPLESHEET_CHECK {
    tag "$samplesheet"

    container "quay.io/wslh-bioinformatics/python@sha256:25b8870fe464a57948723fdc7e8887c003472e2b403997a8da8fdfd1d09b87b1"

    input:
    path samplesheet

    output:
    path '*.csv'       , emit: csv
    path "versions.yml", emit: versions

    script: // This script is bundled with the pipeline, in nf-core/spntypeid/bin/
    """
    check_samplesheet.py \\
        $samplesheet \\
        samplesheet.valid.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
