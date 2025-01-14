#!/usr/bin/env nextflow
/*
========================================================================================
                              Q U A S I F L O W  P I P E L I N E
========================================================================================
              
              A Nextflow pipeline for analysis of NGS-based HIV Drug resitance data
               
----------------------------------------------------------------------------------------
*/

include { runfastQC }   from './modules/runfastQC/main.nf'
include { runMultiQC }  from './modules/runMultiQC/main.nf'
include {} from './modules/'
include {} from './modules/'
include {} from './modules/'
include {} from './modules/'
include {} from './modules/'
include {} from './modules/'
include {} from './modules/'





// Show help emssage
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Header log info
log.info "============================================================"
log.info " AlfredUg/QuasiFlow   ~  version ${params.version}"
log.info "============================================================"
def summary = [:]
summary['Run Name'] = custom_runName ?: workflow.runName
summary['Reads'] = params.reads
summary['Output directory'] = params.outdir
summary['Reporting threshold'] = params.reporting_threshold
summary['Consensus percentage'] = params.consensus_pct
summary['Minimum read quality'] = params.min_read_qual
summary['Length cut off'] = params.length_cutoff
summary['Score cutoff'] = params.score_cutoff
summary['Minimum variant quality'] = params.min_variant_qual
summary['Minimum depth'] = params.min_dp
summary['Minimum allele count'] = params.min_ac
summary['Minimum mutation frequency'] = params.min_freq
summary['Report template'] = params.rmd
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home'] = "$HOME"
summary['Current user'] = "$USER"
summary['Current path'] = "$PWD"
summary['Script dir'] = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(params.email) {
    summary['E-mail Address'] = params.email
}
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { read_pairs_ch_1; read_pairs_ch_2; read_pairs_ch_3 }
 


workflow {

    runfastQC( samples_ch )

    multiqc_ch = runfastQC.out.fastqc_files.collect.()

    runMultiQC( multiqc_ch )



}