process CHECK_EMPTY_NTC {
    label 'process_single'

    container "quay.io/wslh-bioinformatics/pandas@sha256:bf3cb8e5f695cc7c4cf8cc5ab7e7924d1fc4c40dfbe7cb907110e93a7bf6f101"

    input:
    val empty

    output:
    path("Empty_ntcs.tsv"), emit: ntc_samples

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    check_empty_ntc.py \
        -e ${empty}
    """
}
