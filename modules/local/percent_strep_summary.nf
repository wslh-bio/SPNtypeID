process PERCENT_STREP_SUMMARY {
    label 'process_single'

    container "wslh-bioinformatics/pandas@sha256:bf3cb8e5f695cc7c4cf8cc5ab7e7924d1fc4c40dfbe7cb907110e93a7bf6f101"

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
