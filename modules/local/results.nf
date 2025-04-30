process RESULTS {
    label 'process_single'

    container "quay.io/wslh-bioinformatics/pandas@sha256:9ba0a1f5518652ae26501ea464f466dcbb69e43d85250241b308b96406cac458"

    input:
    path("gc_stats_results.tsv")
    path("assembly_stats_results.tsv")
    path("bbduk_results.tsv")
    path("quality_stats.tsv")
    path("coverage_stats.tsv")
    path("quast_results.tsv")
    path("typing_results.tsv")
    path("kraken_ntc_data/*")
    path("kraken_version.yml")

    output:
    path('spntypeid_report.csv')   , emit: result_csv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    results.py
    """
}
