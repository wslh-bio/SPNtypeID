process BBDUK {
    tag "$meta.id"
    label 'process_medium'

    container "quay.io/wslh-bioinformatics/bbtools@sha256:10c3b2d53da7cf6021dfdc684a0119530c4a1aacd9d1a8a2954af2b33440e59e"

    input:
    tuple val(meta), path(reads)
    path(contaminants)

    output:
    tuple val(meta), path('*trimmed*')      , emit: reads
    tuple val(meta), path('*.bbduk.log')    , emit: log
    path("*.trim.txt")                      , emit: bbduk_trim
    path("*.adapter.stats.txt")             , emit: bbduk_adapters
    path "versions.yml"                     , emit: versions
    path('*repaired*')                      , optional: true, emit: repaired_reads
    path('*_singletons.fastq.gz*')          , optional: true, emit: singletons
    path("*.repair.log")                    , optional: true, emit: repaired_log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_type = meta.single_end ? "single" : "paired"
    def repair_in = meta.single_end ? "" : "in1=${reads[0]} in2=${reads[1]}"
    def repair_out = meta.single_end ? "" : "out1=${prefix}_repaired_1.fastq.gz out2=${prefix}_repaired_2.fastq.gz"
    def bbduk_in = meta.single_end ? "in=${reads[0]}" : "in1=${prefix}_repaired_1.fastq.gz in2=${prefix}_repaired_2.fastq.gz"
    def bbduk_out = meta.single_end ? "out=${prefix}_trimmed.fastq.gz" : "out1=${prefix}_trimmed_1.fastq.gz out2=${prefix}_trimmed_2.fastq.gz"
    def contaminants_fa = contaminants ? "ref=$contaminants" : "ref=/bbmap/resources/adapters.fa"
    """
    maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
    read_type=$read_type
    
    if [ \$read_type == "single" ]
    then       
        bbduk.sh \\
            -Xmx\$maxmem \\
            $bbduk_in \\
            $bbduk_out \\
            threads=$task.cpus \\
            $args \\
            $contaminants_fa \\
            stats=${prefix}.adapter.stats.txt \\
            &> ${prefix}.bbduk.log

        grep -E 'Input:|QTrimmed:|Trimmed by overlap:|Total Removed:|Result:' ${prefix}.bbduk.log > ${prefix}.trim.txt
    fi


    if [ \$read_type == "paired" ]
    then
        repair.sh \\
            -Xmx\$maxmem \\
            $repair_in \\
            $repair_out \\
            outs=${prefix}_singletons.fastq.gz \\
            repair \\
            &> ${prefix}.repair.log
        
        bbduk.sh \\
            -Xmx\$maxmem \\
            $bbduk_in \\
            $bbduk_out \\
            threads=$task.cpus \\
            $args \\
            $contaminants_fa \\
            stats=${prefix}.adapter.stats.txt \\
            &> ${prefix}.bbduk.log

        grep -E 'Input:|QTrimmed:|Trimmed by overlap:|Total Removed:|Result:' ${prefix}.bbduk.log > ${prefix}.trim.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh)
    END_VERSIONS
    """
}
