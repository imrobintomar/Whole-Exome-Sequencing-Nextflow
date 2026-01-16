#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.sample_id = null
params.vcf_path = null
params.annovar_path = null
params.hpo_terms = null
params.output_dir = null
params.exomiser_dir = "/path/to/exomiser"
params.exomiser_data_dir = "/path/to/exomiser/data"

include { exomiserPhenotypeAnalysis } from './processes/09_phenotype_analysis.nf'

workflow {
    if (!params.sample_id || !params.vcf_path || !params.annovar_path || !params.hpo_terms) {
        error "Missing required parameters: sample_id, vcf_path, annovar_path, hpo_terms"
    }

    def vcf_file = file(params.vcf_path)
    def annovar_file = file(params.annovar_path)
    def hpo_list = params.hpo_terms.split(',').collect { it.trim() }

    input_ch = Channel.of([params.sample_id, vcf_file, annovar_file, hpo_list])

    exomiserPhenotypeAnalysis(input_ch)
}
