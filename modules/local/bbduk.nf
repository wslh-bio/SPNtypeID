process BBDUK {
    tag "$meta.id"
    label 'process_medium'

    container "wslh-bioinformatics/bbtools@sha256:10c3b2d53da7cf6021dfdc684a0119530c4a1aacd9d1a8a2954af2b33440e59e"

    input:
    tuple val(meta), path(reads)
    path(contaminants)

    output:
    tuple val(meta), path('*.fastq.gz'), emit: reads
    tuple val(meta), path('*.log')     , emit: log
    path("*.trim.txt")                 , emit: bbduk_trim
    path("*.adapter.stats.txt")        , emit: bbduk_adapters
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def raw      = meta.single_end ? "in=${reads[0]}" : "in1=${reads[0]} in2=${reads[1]}"
    def trimmed  = meta.single_end ? "out=${prefix}.fastq.gz" : "out1=${prefix}_trimmed_1.fastq.gz out2=${prefix}_trimmed_2.fastq.gz"
    def contaminants_fa = contaminants ? "ref=$contaminants" : "ref=/bbmap/resources/adapters.fa"
    """
    maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
    bbduk.sh \\
        -Xmx\$maxmem \\
        $raw \\
        $trimmed \\
        threads=$task.cpus \\
        $args \\
        $contaminants_fa \\
        stats=${prefix}.adapter.stats.txt \\
        &> ${prefix}.bbduk.log
    grep -E 'Input:|QTrimmed:|Trimmed by overlap:|Total Removed:|Result:' ${prefix}.bbduk.log > ${prefix}.trim.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh)
    END_VERSIONS
    """
}
