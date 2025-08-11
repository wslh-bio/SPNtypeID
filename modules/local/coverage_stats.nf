process COVERAGE_STATS {
    label 'process_single'

    container "quay.io/wslh-bioinformatics/pandas@sha256:bf3cb8e5f695cc7c4cf8cc5ab7e7924d1fc4c40dfbe7cb907110e93a7bf6f101"

    input:
    path("data*/*")
    val mincoverage

    output:
    path('coverage_stats.tsv'), emit: coverage_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    coverage_stats.py ${mincoverage}
    """
}
