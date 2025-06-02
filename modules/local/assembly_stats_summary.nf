process ASSEMBLY_STATS_SUMMARY {
    label 'process_single'

    container "quay.io/wslh-bioinformatics/spriggan-pandas@sha256:dee64812a5ec4258f9df4b2df47518bd2e7a16cf21fbc354af68874c498b6ce5"

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
