process CREATE_REPORT_NO_NTC {
    label 'process_single'

    container "wslh-bioinformatics/pandas@sha256:bf3cb8e5f695cc7c4cf8cc5ab7e7924d1fc4c40dfbe7cb907110e93a7bf6f101"

    input:
    path asr,           name: "assembly_stats_results.tsv"
    path bb,            name: "bbduk_results.tsv"
    path qs,            name: "quality_stats.tsv"
    path cs,            name: "coverage_stats.tsv"
    path qr,            name: "quast_results.tsv"
    path kv,            name: "kraken_version.yml"
    path psr,           name: "percent_strep_results.tsv"
    path sr,            name: "seroba_results.tsv"
    val run_name_regex
    val split_regex
    val runname


    output:
    path('*_spntypeid_report.csv')   , emit: result_csv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    create_report_no_ntc.py \
        --kraken_version ${kv} \
        --run_name_regex ${run_name_regex} \
        --split_regex ${split_regex} \
        --workflowVersion ${workflow.manifest.version} \
        --workflowRunName ${runname}
    """
}
