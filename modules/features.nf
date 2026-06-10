process FEATURES {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/features", mode: 'copy'

    input:
    tuple val(sample_id), path(fragments)
    path reference
    path regions

    output:
    tuple val(sample_id), path("${sample_id}.features.json"), emit: features

    script:
    """
    fragmentomics features \\
        --fragments ${fragments} \\
        --reference ${reference} \\
        --regions ${regions} \\
        --sample-id ${sample_id} \\
        --k ${params.motif_k} \\
        --out ${sample_id}.features.json
    """
}
