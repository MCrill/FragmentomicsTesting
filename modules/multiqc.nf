process MULTIQC {
    label 'process_low'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path '*'

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc . -n multiqc_report.html
    """
}
