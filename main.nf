#!/usr/bin/env nextflow
/*
 * cfDNA Fragmentomics + NuPEM pipeline
 * FASTQ -> QC/trim -> align -> dedup/filter -> fragments -> feature panel (incl. NuPEM)
 */
nextflow.enable.dsl = 2

include { FASTP        } from './modules/fastp.nf'
include { BWAMEM2      } from './modules/bwamem2.nf'
include { MARKDUP      } from './modules/markdup.nf'
include { FILTER_BAM   } from './modules/filter.nf'
include { EXTRACT_FRAG } from './modules/fragments.nf'
include { FEATURES     } from './modules/features.nf'
include { MULTIQC      } from './modules/multiqc.nf'

workflow {

    // ---- input channels -------------------------------------------------
    // samplesheet: sample_id,fastq_1,fastq_2
    Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, [file(row.fastq_1), file(row.fastq_2)]) }
        .set { reads_ch }

    ref_ch      = file(params.reference,    checkIfExists: true)
    bwa_index   = file(params.bwa_index,    checkIfExists: true)
    regions_ch  = file(params.regions_bed,  checkIfExists: true)
    blacklist   = file(params.blacklist,    checkIfExists: true)

    // ---- workflow -------------------------------------------------------
    FASTP(reads_ch)
    BWAMEM2(FASTP.out.trimmed, bwa_index)
    MARKDUP(BWAMEM2.out.bam)
    FILTER_BAM(MARKDUP.out.bam, blacklist)
    EXTRACT_FRAG(FILTER_BAM.out.bam)
    FEATURES(EXTRACT_FRAG.out.fragments, ref_ch, regions_ch)

    // aggregate QC
    FASTP.out.json
        .mix(MARKDUP.out.metrics)
        .map { it[1] }
        .collect()
        .set { qc_files }
    MULTIQC(qc_files)
}

workflow.onComplete {
    log.info ( workflow.success
        ? "\nDone. Feature tables in ${params.outdir}/features\n"
        : "\nPipeline failed. See ${params.outdir}/pipeline_info\n" )
}
