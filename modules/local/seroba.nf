process SEROBA {
    tag "$meta.id"
    label 'process_medium'

    container "wslh-bioinformatics/seroba@sha256:f805eb8da9e75273589c51a57229d00cbf06ae566ddebd2432f48eee4bcf2614"

    input:
    tuple val(meta), path(reads)

    output:
    path("${meta.id}.pred.csv")                         , emit: seroba_results
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

    # If fastq file ends with _R1_001/_R2_001 or _1/_2 rename those files to end with _R1/_R2 so seroBA does not throw an error
    find . -name "*_R1_001.fastq.gz" -exec bash -c 'mv "\$1" "${prefix}"_R1.fastq.gz' - '{}' +
    find . -name "*_R2_001.fastq.gz" -exec bash -c 'mv "\$1" "${prefix}"_R2.fastq.gz' - '{}' +
    find . -name "*_1.fastq.gz" -exec bash -c 'mv "\$1" "${prefix}"_R1.fastq.gz' - '{}' +
    find . -name "*_2.fastq.gz" -exec bash -c 'mv "\$1" "${prefix}"_R2.fastq.gz' - '{}' +

    # Check to make sure the prefix and file handle match, if not, rename the file with the correct prefix
    for i in "_R1.fastq.gz" "_R2.fastq.gz";
        do
            count_file=`ls -1 *\$i 2>/dev/null | wc -l`
            if [ \$count_file != 0 ]
            then
                f1=`ls *\$i`
                f2="${prefix}""\$i"
                if ! [ "\$f1" = "\$f2" ]
                then
                    mv \$f1 \$f2
                fi
            fi
        done

    seroba runSerotyping $args /seroba*/database ${prefix}_R1.fastq.gz ${prefix}_R2.fastq.gz ${prefix} &> seroba.log
    mv ${prefix}/pred.csv ${prefix}.pred.csv
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
