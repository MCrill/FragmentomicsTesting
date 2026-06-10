process MARKDUP {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/markdup", mode: 'copy', pattern: '*.metrics.txt'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.md.bam"),          emit: bam
    tuple val(sample_id), path("${sample_id}.md.metrics.txt"),  emit: metrics

    script:
    """
    samtools markdup -@ ${task.cpus} \\
        -f ${sample_id}.md.metrics.txt \\
        ${bam} ${sample_id}.md.bam
    samtools index ${sample_id}.md.bam
    """
}
