process CALCULATE_ASSEMBLY_STATS {
    tag"$meta.id"
    container "quay.io/wslh-bioinformatics/spriggan-pandas@sha256:dee64812a5ec4258f9df4b2df47518bd2e7a16cf21fbc354af68874c498b6ce5"

    label 'process_single'

    input:
    tuple val(meta), path(quast_report_tsv)
    path(kraken_results_tsv)
    path NCBI_assembly_stats_file

    output:
    path "*_Assembly_ratio_*"   , emit: assembly_ratio

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in wslh-bio/dryad/bin/
    """
    calculate_assembly_ratio.py \
    -d $NCBI_assembly_stats_file \
    -q $quast_report_tsv \
    -t $kraken_results_tsv
    """
}
