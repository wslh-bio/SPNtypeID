process PERCENT_STREP_SUMMARY {
    label 'process_single'

    container "quay.io/wslh-bioinformatics/pandas@sha256:9ba0a1f5518652ae26501ea464f466dcbb69e43d85250241b308b96406cac458"

    input:
    path("data/*")
    val minpctstrep
    val minpctspn
    val maxpctother

    output:
    path("percent_strep_results.tsv")  , emit: percent_strep_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    percent_strep_summary.py ${minpctstrep} ${minpctspn} ${maxpctother}
    """
}
