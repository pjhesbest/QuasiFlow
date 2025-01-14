process runTrimGalore {

    tag "${sampleID}"

    publishDir "${params.outdir}/adaptors-trimmed-reads", mode: "copy", overwrite: false

    input:
        tuple val(sampleID), 
                path(forward), 
                path(reverse)

    output:
        tuple val(sampleID), 
                path("${sampleID}_val_1.fq.gz"), 
                path("${sampleID}_val_2.fq.gz"),       emit: trimmed_reads_ch

    script:
        """
        trim_galore -q ${params.min_read_qual} \\
                    --basename ${sampleID} \\
                    --paired ${forward} ${reverse} \\
                    -o .

        """
}