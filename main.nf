#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

    /*
    ========================================================================================
                            Q U A S I F L O W  P I P E L I N E
    ========================================================================================
                
                A Nextflow pipeline for analysis of NGS-based HIV Drug resitance data

    ----------------------------------------------------------------------------------------
    */

    /*
        IMPORT MODULES
    */

    include { runfastQC }       from './modules/runfastQC/main.nf'
    include { runMultiQC }      from './modules/runMultiQC/main.nf'
    include { runTrimGalore }   from './modules/runTrimGalore/main.nf'
    include { runHydra }        from './modules/runHydra/main.nf'
    include { runSierralocal }  from './modules/runSierralocal/main.nf'
    include { renderReport }    from './modules/renderReport/main.nf'

    /*
    ······································································································
        REQUIRED ARGUMENTS
    ······································································································
    */

    def helpMessage() {
        log.info"""
        ============================================================
        AlfredUg/QuasiFlow  ~  version ${params.version}
        ============================================================
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run AlfredUg/QuasiFlow --reads <path to fastq files> --outdir <path to output directory>
        
        Mandatory arguments:
            --name                          Name of the run.
            --samplesheet                   Path to input data samplesheet (must be a csv with 4 columns: sampleID,alias,forward_path,reverse_path)

        HyDRA arguments (optional):
            --mutation_db		            Path to mutational database.
            --reporting_threshold	        Minimum mutation frequency percent to report.
            --consensus_pct		            Minimum percentage a base needs to be incorporated into the consensus sequence.
            --min_read_qual	                Minimum quality for a position in a read to be masked.	     
            --length_cutoff	                Reads which fall short of the specified length will be filtered out.
            --score_cutoff		            Reads that have a median or mean quality score (depending on the score type specified) less than the score cutoff value will be filtered out.
            --min_variant_qual              Minimum quality for variant to be considered later on in the pipeline.
            --min_dp                        Minimum required read depth for variant to be considered later on in the pipeline.
            --min_ac                        The minimum required allele count for variant to be considered later on in the pipeline
            --min_freq                      The minimum required frequency for mutation to be considered in drug resistance report.

        Sierralocal arguments (optional):
            --xml                           Path to HIVdb ASI2 XML.
            --apobec-tsv                    Path to tab-delimited (tsv) HIVdb APOBEC DRM file.
            --comments-tsv                  Path to tab-delimited (tsv) HIVdb comments file.
        


        """.stripIndent()
    }

    /*
    ······································································································
        WORKFLOW: main
    ······································································································
    */

    workflow {

        // Create channel from sample sheet
            if (params.samplesheet == null) {
                error "Please provide a samplesheet CSV file with --samplesheet (csv)"
            }

        // Create channel from sample sheet
            if (params.name == null) {
                error "Please provide a runID file with --name (chr)"
            }

    /*
    ······································································································
        CREATION OF CHANNELS
            The section creates the samples_ch and the controls_ch from the samplesheet. First the 
                sample sheet is imported and split by the 'type' column into sample or control, and these
                two sets of samples are directed into seperate workflows.
    ······································································································
    */

        Channel.fromPath(params.samplesheet)
            .splitCsv(header: true, sep: ',')
            .map { row ->
                if (row.sampleID == null || row.forward_path == null || row.reverse_path == null || row.type == null) {
                    error "Missing required column in samplesheet: ${row}"
                }
                tuple( row.sampleID, 
                        file(row.forward_path, checkIfExists: true), 
                        file(row.reverse_path, checkIfExists: true), 
                        )
                }.set { samples_ch }

    /*
    ······································································································
        MAIN WORKFLOW
    ······································································································
    */

        runfastQC( samples_ch )

        multiqc_ch = runfastQC.out.fastqc_files.collect()

        runMultiQC( multiqc_ch )

        runTrimGalore( samples_ch )

        runHydra( runTrimGalore.out.trimmed_reads_ch )
        
        runSierralocal( runHydra.out.cns_sequence_ch )
        
        renderReport( runSierralocal.out.cns_json_ch )        

    }