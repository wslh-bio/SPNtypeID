process TYPING_SUMMARY {
    label 'process_single'

    container "quay.io/wslh-bioinformatics/pandas@sha256:9ba0a1f5518652ae26501ea464f466dcbb69e43d85250241b308b96406cac458"

    input:
    path("data/*")
    val minpctstrep
    val minpctspn
    val maxpctother

    output:
    path("typing_results.tsv")  , emit: typing_summary_results
    path("kraken_results.tsv")
    path("seroba_results.tsv")

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    typing_summary.py ${minpctstrep} ${minpctspn} ${maxpctother}
    """
}
