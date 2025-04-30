process QUALITY_STATS {
    label 'process_single'

    container "quay.io/wslh-bioinformatics/pandas@sha256:9ba0a1f5518652ae26501ea464f466dcbb69e43d85250241b308b96406cac458"

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
