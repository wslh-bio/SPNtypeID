process KRAKEN_SUMMARY {
    label 'process_single'

    container "quay.io/wslh-bioinformatics/pandas@sha256:9ba0a1f5518652ae26501ea464f466dcbb69e43d85250241b308b96406cac458"

    input:
    path("data/*")

    output:
    path("kraken_results.tsv"), emit: kraken_tsv

    script:
    """
    kraken_summary.py
    """
}
