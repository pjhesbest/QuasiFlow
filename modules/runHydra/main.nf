process runHydra{

    tag "$sampleID"
    
    publishDir params.outdir, mode: 'copy'

    input:
        tuple val(sampleID), 
                    path("${sampleID}_val_1.fq.gz"), 
                    path("${sampleID}_val_2.fq.gz")
    
    output:
        tuple val(sampleID), path('consensus_*.fasta'),         emit: cns_sequence_ch
        tuple val(sampleID), path('dr_report_*.csv'),           emit: dr_report_ch
        tuple val(sampleID), path('mutation_report_*.aavf'),    emit: dr_report_ch_2
        tuple val(sampleID), path('filtered_*.fastq'),          emit: filtered_ch

    script:
        """
        # Run quasitools
            quasitools hydra \\
                ${sampleID}_val_1.fq.gz ${sampleID}_val_2.fq.gz \\
                -o . \\
                --generate_consensus \\
                --reporting_threshold ${params.reporting_threshold} \\
                --consensus_pct ${params.consensus_pct} \\
                --length_cutoff ${params.length_cutoff} \\
                --score_cutoff ${params.score_cutoff} \\
                --min_variant_qual ${params.min_variant_qual} \\
                --min_dp ${params.min_dp} \\
                --min_ac ${params.min_ac} \\
                --min_freq ${params.min_freq}
        
        # rename to ensure results are unqiue
            mv consensus.fasta consensus_${sampleID}.fasta 
            mv dr_report.csv dr_report_${sampleID}.csv
            mv mutation_report.aavf mutation_report_${sampleID}.aavf
            mv filtered.fastq filtered_${sampleID}.fastq

        """
}