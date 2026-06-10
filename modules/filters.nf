process FILTER_BAM {
    tag "$sample_id"
    label 'process_low'

    input:
    tuple val(sample_id), path(bam)
    path blacklist

    output:
    tuple val(sample_id), path("${sample_id}.filt.bam"), emit: bam

    script:
    // keep properly-paired, primary, non-dup, MAPQ>=30 reads on autosomes,
    // then drop reads overlapping the ENCODE blacklist.
    """
    samtools view -@ ${task.cpus} -b \\
        -f 0x2 -F 0xF04 -q ${params.min_mapq} \\
        ${bam} \$(seq 1 22 | sed 's/^/chr/') \\
    | bedtools intersect -v -abam stdin -b ${blacklist} \\
    > ${sample_id}.filt.bam
    samtools index ${sample_id}.filt.bam
    """
}
