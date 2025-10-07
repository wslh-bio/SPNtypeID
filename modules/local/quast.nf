process QUAST {
    tag "$meta.id"
    label 'process_medium'

    container "wslh-bioinformatics/quast@sha256:84753d0c00487d7f5e99265f4912efeb0823c9b86731e89cadd6384f913c1c05"

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("${meta.id}.transposed.quast.report.tsv") , emit: transposed_report
    path("${meta.id}.quast.report.tsv")             , emit: result
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    quast.py ${contigs} -o .
    mv report.tsv ${prefix}.quast.report.tsv
    mv transposed_report.tsv ${prefix}.transposed.quast.report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast.py --version 2>&1 | sed 's/^.*QUAST v//; s/ .*\$//')
    END_VERSIONS
    """
}
