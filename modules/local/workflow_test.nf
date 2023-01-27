process WORKFLOW_TEST {
    tag "workflow_validation"
    label 'process_low'

    container "quay.io/wslh-bioinformatics/pandas@sha256:bf3cb8e5f695cc7c4cf8cc5ab7e7924d1fc4c40dfbe7cb907110e93a7bf6f101"

    input:
    path(valid_data)
    path(test_data)

    output:
    path "validation.log"

    when:
    workflow_test = true

    script:
    """
    workflow_validation.py $valid_data $test_data | tee validation.log
    """
}
