process CREATE_REPORT {
    label 'process_single'

    container "quay.io/wslh-bioinformatics/pandas@sha256:bf3cb8e5f695cc7c4cf8cc5ab7e7924d1fc4c40dfbe7cb907110e93a7bf6f101"

    input:
    path results_compiled
    val empty_ntc
    val runname

    output:
    path('*_spntypeid_report.csv')   , emit: result_csv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    create_report.py \
        --result_files ${results_compiled} \
        --workflowVersion ${workflow.manifest.version} \
        --workflowRunName ${runname} \
        --empty_ntc_list ${empty_ntc}
    """
}
