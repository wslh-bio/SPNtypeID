process CALCULATE_ASSEMBLY_STATS {
    tag "$meta.id"
    label 'process_single'

    container "quay.io/wslh-bioinformatics/pandas@sha256:bf3cb8e5f695cc7c4cf8cc5ab7e7924d1fc4c40dfbe7cb907110e93a7bf6f101"

    input:
    tuple val(meta), path(quast_report_tsv)
    path NCBI_assembly_stats_file


    output:
    path "*_Assembly_ratio_*"   , emit: assembly_ratio

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in wslh-bio/dryad/bin/
    """
    calculate_assembly_ratio.py \
        -q $quast_report_tsv \
        -d $NCBI_assembly_stats_file

    """
}
