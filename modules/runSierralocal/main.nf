process runSierralocal {

    conda params.sierra_local_env

    tag "$sampleID"

    publishDir params.outdir, mode: 'copy'

    input:
        tuple val(sampleID), path(cns_seq)

    output:
        tuple val(sampleID), path('consensus*.json'),   emit: cns_json_ch

    script:

        """
        sierralocal ${cns_seq} -o consensus_${sampleID}.json -xml ${params.xml}

        """
}