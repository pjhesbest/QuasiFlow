process runMultiQC{
    
    publishDir "${params.outdir}/qcresults-raw-reads", mode: "copy", overwrite: false

    input:
        path(multiqc_ch)

    output:
        path('raw_reads_multiqc_report.html')

    script:
        """

        multiqc .

        mv multiqc_report.html raw_reads_multiqc_report.html

        """
}
