process QUAST_SUMMARY {
    label 'process_single'

    container "quay.io/wslh-bioinformatics/pandas@sha256:bf3cb8e5f695cc7c4cf8cc5ab7e7924d1fc4c40dfbe7cb907110e93a7bf6f101"

    input:
    path("data*/*")
    val maxcontigs

    output:
    path("quast_results.tsv"), emit: quast_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    quast_summary.py ${maxcontigs}
    """
}
