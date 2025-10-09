process CREATE_REPORT {
    label 'process_single'

    container "quay.io/wslh-bioinformatics/pandas@sha256:bf3cb8e5f695cc7c4cf8cc5ab7e7924d1fc4c40dfbe7cb907110e93a7bf6f101"

    input:
    path asr,           name: "assembly_stats_results.tsv"
    path bb,            name: "bbduk_results.tsv"
    path qs,            name: "quality_stats.tsv"
    path cs,            name: "coverage_stats.tsv"
    path qr,            name: "quast_results.tsv"
    path kntc,          name: "kraken_ntc_data/*"
    path kv,            name: "kraken_version.yml"
    path psr,           name: "percent_strep_results.tsv"
    path sr,            name: "seroba_results.tsv"
    val en
    val runname

    output:
    path('*_spntypeid_report.csv')   , emit: result_csv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    create_report.py \
        --kraken_ntc_data ${kntc} \
        --kraken_version ${kv} \
        --workflowVersion ${workflow.manifest.version} \
        --workflowRunName ${runname} \
        --empty_ntc_list ${en}
    """
}
