process BWAMEM2 {
    tag "$sample_id"
    label 'process_high'

    input:
    tuple val(sample_id), path(reads)
    path index

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), emit: bam

    script:
    // index is the directory/prefix produced by `bwa-mem2 index`
    def prefix = index.find { it.name.endsWith('.bwt.2bit.64') }?.simpleName ?: index[0].simpleName
    """
    bwa-mem2 mem -t ${task.cpus} \\
        -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" \\
        ${prefix} ${reads[0]} ${reads[1]} \\
    | samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam -
    samtools index ${sample_id}.sorted.bam
    """
}
