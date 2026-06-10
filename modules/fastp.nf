process FASTP {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/fastp", mode: 'copy', pattern: '*.{json,html}'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_{1,2}.trim.fastq.gz"), emit: trimmed
    tuple val(sample_id), path("${sample_id}.fastp.json"),          emit: json
    path "${sample_id}.fastp.html",                                 emit: html

    script:
    """
    fastp \\
        -i ${reads[0]} -I ${reads[1]} \\
        -o ${sample_id}_1.trim.fastq.gz -O ${sample_id}_2.trim.fastq.gz \\
        --detect_adapter_for_pe \\
        --qualified_quality_phred 20 \\
        --length_required 30 \\
        --thread ${task.cpus} \\
        --json ${sample_id}.fastp.json \\
        --html ${sample_id}.fastp.html
    """
}
