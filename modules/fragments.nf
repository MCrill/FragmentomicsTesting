process EXTRACT_FRAG {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.outdir}/fragments", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.fragments.bed.gz"), emit: fragments
    path "${sample_id}.fragments.bed.gz.tbi"

    script:
    """
    fragmentomics extract \\
        --bam ${bam} \\
        --out ${sample_id}.fragments.bed \\
        --min-mapq ${params.min_mapq} \\
        --min-length ${params.min_frag_len} \\
        --max-length ${params.max_frag_len}
    sort -k1,1 -k2,2n ${sample_id}.fragments.bed | bgzip > ${sample_id}.fragments.bed.gz
    tabix -p bed ${sample_id}.fragments.bed.gz
    """
}
