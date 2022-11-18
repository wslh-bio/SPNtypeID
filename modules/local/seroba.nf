process SEROBA {
    tag "$meta.id"
    label 'process_medium'

    container "quay.io/wslh-bioinformatics/seroba@sha256:1846871901c1dfe25c43028ed0ec326bb0747d5eb11a5d907b3fb35197c420fa"

    input:
    tuple val(meta), path(reads)

    output:
    path("${meta.id}.pred.tsv")                         , emit: seroba_results
    path("${meta.id}_detailed_serogroup_info.txt")      , optional: true, emit: seroba_detailed_results
    path("seroba.log")                                  , emit: log
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def memory = task.memory.toGiga()
    """
    export MPLCONFIGDIR=\$(pwd)

    ln -s ${reads[0]} ${prefix}_R1.fastq.gz
    ln -s ${reads[1]} ${prefix}_R2.fastq.gz 

    seroba runSerotyping $args /seroba*/database ${prefix}_R1.fastq.gz ${prefix}_R2.fastq.gz ${prefix} &> seroba.log
    mv ${prefix}/pred.tsv ${prefix}.pred.tsv
    if [ -f detailed_serogroup_info.txt ]; then
        mv detailed_serogroup_info.txt ${prefix}_detailed_serogroup_info.txt
    fi
    find -name  seroba.log -size 0 -exec rm {} +

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seroba: \$(echo \$(seroba version 2>&1))
    END_VERSIONS
    """
}
