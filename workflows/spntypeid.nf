/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//


include { BBDUK                         } from '../modules/local/bbduk'
include { BBDUK_SUMMARY                 } from '../modules/local/bbduk_summary'
include { FASTQC                        } from '../modules/local/fastqc'
include { FASTQC_SUMMARY                } from '../modules/local/fastqc_summary'
include { SHOVILL                       } from '../modules/local/shovill'
include { SAMTOOLS                      } from '../modules/local/samtools'
include { COVERAGE_STATS                } from '../modules/local/coverage_stats'
include { QUAST                         } from '../modules/local/quast'
include { QUAST_SUMMARY                 } from '../modules/local/quast_summary'
include { BIOAWK                        } from '../modules/local/bioawk'
include { QUALITY_STATS                 } from '../modules/local/quality_stats'
include { KRAKEN as KRAKEN_SAMPLE       } from '../modules/local/kraken'
include { KRAKEN as KRAKEN_NTC          } from '../modules/local/kraken'
include { SEROBA                        } from '../modules/local/seroba'
include { TYPING_SUMMARY                } from '../modules/local/typing_summary'
include { RESULTS                       } from '../modules/local/results'
include { WORKFLOW_TEST                 } from '../modules/local/workflow_test'
include { MULTIQC                       } from '../modules/local/multiqc'
include { CUSTOM_DUMPSOFTWAREVERSIONS   } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow SPNTYPEID {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )

    INPUT_CHECK.out.reads
        .map { meta, file ->
        [meta, file, file[0].countFastq()]
        }
        .branch{ meta, file, count ->
            pass: count > 0
            fail: count == 0
        }
        .set{ ch_count }

    ch_count.pass
        .map { meta, file, count ->
            [meta, file]
            }
        .set{ ch_filtered }

    ch_count.fail
        .map{ meta, file, count -> 
            [meta.id]
        }
        .flatten()
        .set{ch_failed}
    ch_failed
        .collectFile(
            storeDir: "${params.outdir}/rejected_samples",
            name: 'Empty_samples.csv',
            newLine: true
        )

    ch_filtered
        .branch {
            ntc: !!(it[0]['id'] =~ params.ntc_regex)
            sample: true
        }
        .set{ ch_input_reads }

    //
    // MODULE: BBDUK
    //
    BBDUK (
        ch_input_reads.sample,
        params.contaminants
    )
    ch_versions = ch_versions.mix(BBDUK.out.versions.first())

    //
    // MODULE: BBDUK_SUMMARY
    //
    BBDUK_SUMMARY (
        BBDUK.out.bbduk_trim.collect()
    )

    //
    // MODULE: FASTQC
    //
    FASTQC (
        BBDUK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: FASTQC_SUMMARY
    //
    FASTQC_SUMMARY (
        FASTQC.out.zip.collect{it[1]}
    )

    //
    // MODULE: SHOVILL
    //
    SHOVILL (
        BBDUK.out.reads
    )
    ch_versions = ch_versions.mix(SHOVILL.out.versions.first())

    //
    // MODULE: SAMTOOLS
    //
    SAMTOOLS (
        SHOVILL.out.sam_files
    )
    ch_versions = ch_versions.mix(SAMTOOLS.out.versions.first())

    //
    // MODULE: COVERAGE_STATS
    //
    COVERAGE_STATS (
        SAMTOOLS.out.cov_files.collect()
    )

    //
    // MODULE: QUAST
    //
    QUAST (
        SHOVILL.out.contigs
    )
    ch_versions = ch_versions.mix(QUAST.out.versions.first())

    //
    // MODULE: QUAST_SUMMARY
    //
    QUAST_SUMMARY (
        QUAST.out.transposed_report.collect()
    )

    //
    // MODULE: BIOAWK
    //
    BIOAWK (
        ch_input_reads.sample
    )
    ch_versions = ch_versions.mix(BIOAWK.out.versions.first())

    //
    // MODULE: QUALITY_STATS
    //
    QUALITY_STATS (
        BIOAWK.out.qual_results.collect()
    )

    //
    // MODULE: KRAKEN_SAMPLE
    //
    KRAKEN_SAMPLE (
        ch_input_reads.sample
    )
    ch_versions = ch_versions.mix(KRAKEN_SAMPLE.out.versions.first())

    //
    // MODULE: KRAKEN_NTC
    //
    KRAKEN_NTC (
        ch_input_reads.ntc
    )

    //
    // MODULE: SEROBA
    //
    SEROBA (
        ch_input_reads.sample
    )
    ch_versions = ch_versions.mix(SEROBA.out.versions.first())

    //
    // MODULE: TYPING_SUMMARY
    //
    TYPING_SUMMARY (
        KRAKEN_SAMPLE.out.kraken_results.mix(SEROBA.out.seroba_results).collect()
    )

    //
    // MODULE: RESULTS
    //
    RESULTS (
        BBDUK_SUMMARY.out.bbduk_tsv,
        QUALITY_STATS.out.quality_tsv,
        COVERAGE_STATS.out.coverage_tsv,
        QUAST_SUMMARY.out.quast_tsv,
        TYPING_SUMMARY.out.typing_summary_results,
        KRAKEN_NTC.out.kraken_results.collect().ifEmpty([]),
        KRAKEN_SAMPLE.out.versions.first()
    )

    //
    // MODULE: WORKFLOW_TEST
    //
    ch_valid_dataset = Channel.fromPath("$projectDir/test-dataset/validation/spntypeid_report_valid.csv", checkIfExists: true)
    WORKFLOW_TEST (
        ch_valid_dataset.collect(),
        RESULTS.out.result_csv
    )


    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowSpntypeid.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowSpntypeid.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BBDUK.out.bbduk_adapters.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BBDUK.out.bbduk_trim.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS.out.stats_multiqc.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(KRAKEN_SAMPLE.out.kraken_results.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(KRAKEN_NTC.out.kraken_results.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.result.collect().ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.collect().ifEmpty([]),
        ch_multiqc_custom_config.collect().ifEmpty([]),
        ch_multiqc_logo.collect().ifEmpty([])
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.adaptivecard(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/