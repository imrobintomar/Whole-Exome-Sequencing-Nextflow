#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ===========================================
//
// Author: Robin Tomar
// Email: Aiimsgenomics@gmail.com
// GitHub: https://github.com/imrobintomar/Whole-Exome-Sequencing-Nextflow.git
// ===========================================

// Script parameters
params.input_dir = "/media/drprabudh/m3/PRJNA855946/FASTQ"
params.reference = "/media/drprabudh/m1/hg38/hg38.fa"
params.known_sites = [
    "/media/drprabudh/m1/vcf_file/Homo_sapiens_assembly38.known_indels.vcf.gz",
    "/media/drprabudh/m1/vcf_file/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
]
params.output_dir = "/media/drprabudh/m3/PRJNA855946"
params.intervals = "/media/drprabudh/m1/Downloads/WES_Final_Agilent_V8_Hg38_1.bed"  // Default exome target regions
params.help = false

// Include processes
include { fastpQC } from './processes/01_qc.nf'
include { bwaMem } from './processes/02_alignment.nf'
include { sortSam } from './processes/03_sorting.nf'
include { flagstat } from './processes/03_sorting.nf'
include { markDuplicates } from './processes/04_dedup.nf'
include { sortSamPostDedup } from './processes/04_dedup.nf'
include { baseRecalibrator } from './processes/05_bqsr.nf'
include { applyBQSR } from './processes/05_bqsr.nf'
include { haplotypeCaller } from './processes/06_variantcalling.nf'
include { annovarAnnotate } from './processes/07_annotation.nf'
include { addUniqueID } from './processes/08_filtering.nf'

workflow {
    // Additional runtime validation check
    def requiredAuthor = 'Robin Tomar'
    def requiredGithub = 'https://github.com/imrobintomar/Whole-Exome-Sequencing-Nextflow.git'

    if (workflow.manifest.author != requiredAuthor || workflow.manifest.homePage != requiredGithub) {
        error """

        ═══════════════════════════════════════════════════════════════
                        PIPELINE EXECUTION BLOCKED
        ═══════════════════════════════════════════════════════════════

        ERROR: Author attribution has been tampered with!

        This pipeline will NOT run with modified author information.

        Required:
          Author:  ${requiredAuthor}
          GitHub:  ${requiredGithub}

        Found:
          Author:  ${workflow.manifest.author ?: 'MISSING'}
          GitHub:  ${workflow.manifest.homePage ?: 'MISSING'}

        Contact: Aiimsgenomics@gmail.com
        ═══════════════════════════════════════════════════════════════
        """
    }

    if (params.help) {
        println """
        ==========================================
        Whole Exome Sequencing Pipeline
        ==========================================
        Author:  Robin Tomar
        Email:   Aiimsgenomics@gmail.com
        GitHub:  https://github.com/imrobintomar/Whole-Exome-Sequencing-Nextflow.git

        Usage:
            nextflow run main.nf --input_dir <path> --reference <path> --output_dir <path>

        Options:
            --input_dir         Input FASTQ directory (default: ${params.input_dir})
            --reference         Reference genome FASTA (default: ${params.reference})
            --output_dir        Output directory (default: ${params.output_dir})
            --intervals         Exome target regions BED file (default: ${params.intervals})
            --help              Show this help message

        Resume Failed Pipeline:
            nextflow run main.nf -resume

        Clean Restart:
            nextflow run main.nf -with-trace
        """
        exit 0
    }

    // Create output directories
    file(params.output_dir).mkdirs()
    file("${params.output_dir}/filtered_fastp").mkdirs()
    file("${params.output_dir}/Mapsam").mkdirs()
    file("${params.output_dir}/Germline_VCF").mkdirs()
    file("${params.output_dir}/logs").mkdirs()

    // ═══════════════════════════════════════════════════════════════
    // CRITICAL: Validate all required files exist before execution
    // ═══════════════════════════════════════════════════════════════

    // 1. Validate reference FASTA and indices
    if (!file(params.reference).exists()) {
        error """
        ═══════════════════════════════════════════════════════════════
        VALIDATION ERROR: Reference FASTA not found!
        ═══════════════════════════════════════════════════════════════

        Missing file: ${params.reference}

        Please ensure the reference genome FASTA file exists.
        ═══════════════════════════════════════════════════════════════
        """
    }

    if (!file("${params.reference}.fai").exists()) {
        error """
        ═══════════════════════════════════════════════════════════════
        VALIDATION ERROR: Reference index (.fai) not found!
        ═══════════════════════════════════════════════════════════════

        Missing file: ${params.reference}.fai

        Please create the index by running:
        samtools faidx ${params.reference}
        ═══════════════════════════════════════════════════════════════
        """
    }

    def dict_file = params.reference_dict ?: params.reference.replaceAll(/\.fa(sta)?$/, '.dict')
    if (!file(dict_file).exists()) {
        error """
        ═══════════════════════════════════════════════════════════════
        VALIDATION ERROR: Reference dictionary (.dict) not found!
        ═══════════════════════════════════════════════════════════════

        Missing file: ${dict_file}

        Please create the dictionary by running:
        gatk CreateSequenceDictionary -R ${params.reference}
        ═══════════════════════════════════════════════════════════════
        """
    }

    // 2. Validate known sites VCF files and indices
    params.known_sites.each { vcf ->
        if (!file(vcf).exists()) {
            error """
            ═══════════════════════════════════════════════════════════════
            VALIDATION ERROR: Known sites VCF not found!
            ═══════════════════════════════════════════════════════════════

            Missing file: ${vcf}

            Please ensure all known sites VCF files exist.
            ═══════════════════════════════════════════════════════════════
            """
        }
        if (!file("${vcf}.tbi").exists()) {
            error """
            ═══════════════════════════════════════════════════════════════
            VALIDATION ERROR: Known sites VCF index (.tbi) not found!
            ═══════════════════════════════════════════════════════════════

            Missing file: ${vcf}.tbi

            Please create the index by running:
            tabix -p vcf ${vcf}
            ═══════════════════════════════════════════════════════════════
            """
        }
    }

    // 3. Validate intervals file if provided (recommended for exome)
    if (params.intervals && !file(params.intervals).exists()) {
        error """
        ═══════════════════════════════════════════════════════════════
        VALIDATION ERROR: Intervals file not found!
        ═══════════════════════════════════════════════════════════════

        Missing file: ${params.intervals}

        Please ensure the intervals/BED file exists, or set:
        --intervals null

        to analyze the whole genome (not recommended for exome).
        ═══════════════════════════════════════════════════════════════
        """
    }

    // 4. Warn if intervals not set for exome analysis
    if (!params.intervals) {
        log.warn """
        ═══════════════════════════════════════════════════════════════
        WARNING: No intervals file specified!
        ═══════════════════════════════════════════════════════════════

        For EXOME sequencing, you should specify target regions:
        --intervals /path/to/exome_targets.bed

        Without intervals:
        - Analysis will be 60x SLOWER
        - False positive rate will be HIGHER
        - Resource usage will be EXCESSIVE

        Only proceed if analyzing whole genome data.
        ═══════════════════════════════════════════════════════════════
        """
    }

    // Log pipeline start
    log.info """
    ╔════════════════════════════════════════╗
    ║   Whole Exome Sequencing Pipeline      ║
    ║                                        ║
    ╚════════════════════════════════════════╝
    ==========================================
        Author:  Robin Tomar
        Email:   Aiimsgenomics@gmail.com
    ==========================================


    Input Directory:     ${params.input_dir}
    Output Directory:    ${params.output_dir}
    Reference Genome:    ${params.reference}
    Intervals File:      ${params.intervals ?: 'NOT SET (whole genome mode)'}

    Pipeline started at: ${new Date()}
    """

    // Check if filtered files already exist
    def filtered_dir = "${params.output_dir}/filtered_fastp"
    def filtered_files = file("${filtered_dir}/*_filtered.fastq.gz")
    def filtered_exists = filtered_files.size() > 0

    if (filtered_exists) {
        log.info "✓ Found existing filtered files in ${filtered_dir}, skipping fastpQC step"
        // Read already-filtered paired-end FASTQ files
        fastp_out = Channel.fromFilePairs("${filtered_dir}/*_{1,2}_filtered.fastq.gz", size: 2)
            .ifEmpty { error "No filtered FASTQ files found in ${filtered_dir}" }
            .map { sample_id, files ->
                // Extract base sample ID (remove _1_filtered or _2_filtered suffix)
                def base_id = sample_id.replaceAll(/_[12]$/, '')
                tuple(base_id, files[0], files[1])
            }
    } else {
        log.info "Running fastpQC step on raw FASTQ files"
        // Read paired-end FASTQ files
        read_pairs = Channel.fromFilePairs("${params.input_dir}/*_{1,2}.fastq.gz", size: 2)
            .ifEmpty { error "No FASTQ files found in ${params.input_dir}" }

        // QC with Fastp
        fastp_out = fastpQC(read_pairs)
    }

    // BWA alignment
    bwa_out = bwaMem(fastp_out)

    // Sort SAM to BAM
    bam_sorted = sortSam(bwa_out)

    // Flagstat
    flagstat(bam_sorted)

    // Mark duplicates
    marked_bam = markDuplicates(bam_sorted)
    marked_sorted = sortSamPostDedup(marked_bam)

    // BQSR
    recal_table = baseRecalibrator(marked_sorted)
    final_bam = applyBQSR(marked_sorted, recal_table)

    // Variant calling
    vcf_raw = haplotypeCaller(final_bam)

    // ANNOVAR annotation (single source of truth)
    annovar_txt = annovarAnnotate(vcf_raw)

    // Add UniqueID to ANNOVAR output
    final_output = addUniqueID(annovar_txt)

    // Pipeline completion message
    final_output.subscribe {
        log.info "✓ Pipeline completed successfully at ${new Date()}"
    }
}