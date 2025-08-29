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
    path en,            name: "Empty_ntcs.tsv"
    val ntc_read_limit
    val ntc_spn_read_limit
    val run_name_regex
    val split_regex
    val min_assembly_length
    val max_assembly_length

    output:
    path('*_spntypeid_report.csv')   , emit: result_csv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    create_report.py \
        --kraken_ntc_data ${kntc} \
        --kraken_version ${kv} \
        --ntc_read_limit ${ntc_read_limit} \
        --ntc_spn_read_limit ${ntc_spn_read_limit} \
        --run_name_regex ${run_name_regex} \
        --split_regex ${split_regex} \
        --min_assembly_length ${min_assembly_length} \
        --max_assembly_length ${max_assembly_length} \
        --workflowVersion ${workflow.manifest.version} \
        --workflowRunName ${workflow.runName}
        --empty_ntc_file ${en}
    """
}
