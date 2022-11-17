process BIOAWK {
    tag "$meta.id"
    label 'process_single'

    container "quay.io/wslh-bioinformatics/bioawk@sha256:9e45d66b690722172f1c0e3396cc93373ade3bbe75c60cae338a0eb3ac2f44c9"

    input:
    tuple val(meta), path(reads)

    output:
    path("*.qual.tsv")          , emit: qual_results
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat ${reads[0]} ${reads[1]} > ${prefix}_qc.fastq.gz
    bioawk -c fastx '{print ">"\$name; print meanqual(\$qual)}' ${prefix}_qc.fastq.gz > ${prefix}.bioawk.tsv
    awk 'NR % 2 == 0' ${prefix}.bioawk.tsv > ${prefix}.qual.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioawk: \$(echo \$(bioawk --version 2>&1))
    END_VERSIONS
    """
}
