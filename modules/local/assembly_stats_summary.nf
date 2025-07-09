process ASSEMBLY_STATS_SUMMARY {
    label 'process_single'

    container "quay.io/wslh-bioinformatics/pandas@sha256:bf3cb8e5f695cc7c4cf8cc5ab7e7924d1fc4c40dfbe7cb907110e93a7bf6f101"

    input:
    path("data/*")

    output:
    path("assembly_stats_results_summary.tsv"), emit: assembly_stats_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    assembly_stats_summary.py
    """
}
