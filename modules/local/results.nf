process RESULTS {
    label 'process_single'

    container "quay.io/wslh-bioinformatics/pandas@sha256:9ba0a1f5518652ae26501ea464f466dcbb69e43d85250241b308b96406cac458"

    input:
    path gc,            name: "gc_stats_results.tsv"
    path asr,           name: "assembly_stats_results.tsv"
    path bb,            name: "bbduk_results.tsv"
    path qs,            name: "quality_stats.tsv"
    path cs,            name: "coverage_stats.tsv"
    path qr,            name: "quast_results.tsv"
    path kntc,          name: "kraken_ntc_data/*"
    path kv,            name: "kraken_version.yml"
    path psr,           name: "percent_strep_results.tsv"
    path sr,            name: "seroba_results.tsv"
    val ntc_read_limit
    val ntc_spn_read_limit
    val run_name_regex
    val split_regex

    output:
    path('*_spntypeid_report.csv')   , emit: result_csv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    results.py \
        --gc_stats ${gc} \
        --assembly_stats ${asr} \
        --bbduk_stats ${bb} \
        --quality_stats ${qs} \
        --coverage_stats ${cs} \
        --quast_results ${qr} \
        --kraken_ntc_data ${kntc} \
        --kraken_version ${kv} \
        --ntc_read_limit ${ntc_read_limit} \
        --ntc_spn_read_limit ${ntc_spn_read_limit} \
        --workflowVersion ${workflow.manifest.version} \
        --workflowRunName ${workflow.runName} \
        --percentStrepResults ${psr} \
        --serobaResults ${sr}
    """
}