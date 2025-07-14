process BBDUK_SUMMARY {
    label 'process_single'

    container "wslh-bioinformatics/pandas@sha256:bf3cb8e5f695cc7c4cf8cc5ab7e7924d1fc4c40dfbe7cb907110e93a7bf6f101"

    input:
    path("data*/*")

    output:
    path("bbduk_results.tsv"), emit: bbduk_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bbduk_summary.py 
    """
}
