process WORKFLOW_TEST {
    tag "workflow_validation"
    label 'process_low'

    container "quay.io/wslh-bioinformatics/pandas@sha256:9ba0a1f5518652ae26501ea464f466dcbb69e43d85250241b308b96406cac458"

    input:
    path(valid_data)
    path(test_data)

    output:
    path "validation.log"

    when:
    workflow_test = true

    script:
    """
    workflow_validation.py $valid_data $test_data | tee > validation.log
    """
}
