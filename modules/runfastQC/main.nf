process runfastQC {

        tag "${sampleID}" 

        publishDir "${params.outdir}/qcresults-raw-reads", mode: "copy"

        input:
                tuple val(sampleID), 
                        path(forward), 
                        path(reverse)

        output:
        path("${sampleID}_fastqc/*.zip"), emit: fastqc_files

        script:
                """
                mkdir ${sampleID}_fastqc

                fastqc --outdir ${sampleID}_fastqc \\
                        ${forward} \\
                        ${reverse}

        """
}