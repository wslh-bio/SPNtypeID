params.workflow_test = false
process WORKFLOW_TEST {
    tag "workflow_validation"
    label 'process_low'

    container "quay.io/wslh-bioinformatics/pandas@sha256:bf3cb8e5f695cc7c4cf8cc5ab7e7924d1fc4c40dfbe7cb907110e93a7bf6f101"

    input:
    path(valid_data)
    path(test_data)
    val (sample_1_z_avg)
    val (sample_1_ratio_avg)
    val (sample_2_z_avg)
    val (sample_2_ratio_avg)
    val (sample_3_z_avg)
    val (sample_3_ratio_avg)

    output:
    path "validation.log"

    when:
    params.workflow_test == true

    script:
    """
    workflow_validation.py $valid_data $test_data --sample_1_z_avg $sample_1_z_avg --sample_1_ratio_avg $sample_1_ratio_avg --sample_2_z_avg $sample_2_z_avg --sample_2_ratio_avg $sample_2_ratio_avg --sample_3_z_avg $sample_3_z_avg --sample_3_ratio_avg $sample_3_ratio_avg | tee validation.log
    """
}
