process QUALITY_STATS {
    label 'process_single'

    container "wslh-bioinformatics/pandas@sha256:bf3cb8e5f695cc7c4cf8cc5ab7e7924d1fc4c40dfbe7cb907110e93a7bf6f101"

    input:
    path("data*/*")
    val minavgreadq

    output:
    path('quality_stats.tsv')   , emit: quality_tsv


    when:
    task.ext.when == null || task.ext.when

    script:
    """
    quality_stats.py ${minavgreadq}
    """
}
