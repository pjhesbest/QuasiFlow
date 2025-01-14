process renderReport{

    tag "$sampleID"

    publishDir params.outdir, mode: 'copy'
    
    input:
        tuple val(sampleID), path(cns_json)

    output:
        tuple val(sampleID), path('hivdr_*.html'), path('hivdr_*.pdf'), emit: final_report_ch

    script:

        // path rmd from params.rmd
        // path rmd_static from params.rmd_static 

        """
        Rscript -e 'rmarkdown::render("/home/phesketh/Documents/GitHub/QuasiFlow/assets/hivdr.Rmd", 
            params=list(
                mutation_comments="${params.mutation_db_comments}", 
                dr_report_hivdb="${params.outdir}/${cns_json}",
                mutational_threshold=${params.min_freq},
                minimum_read_depth=${params.min_dp},
                minimum_percentage_cons=${params.consensus_pct}), 
                output_file="hivdr_${sampleID}.html", output_dir = getwd())'

        Rscript -e 'rmarkdown::render("/home/phesketh/Documents/GitHub/QuasiFlow/assets/hivdr_static.Rmd", 
            params=list(
                mutation_comments="${params.mutation_db_comments}",
                dr_report_hivdb="${params.outdir}/${cns_json}",
                mutational_threshold=${params.min_freq},
                minimum_read_depth=${params.min_dp},
                minimum_percentage_cons=${params.consensus_pct}), 
                output_file="hivdr_${sampleID}.pdf", output_dir = getwd())'

        """
}