process baseRecalibrator {

    tag "$sample_id"

    input:
        tuple val(sample_id), path(bam_file)

    output:
        tuple val(sample_id), path("${sample_id}_recall.table")

    script:
        def known_sites_str = params.known_sites
                                  .collect { "--known-sites ${it}" }
                                  .join(" ")

        """
        gatk BaseRecalibrator \
            -R ${params.reference} \
            -I ${bam_file} \
            ${known_sites_str} \
            -O ${sample_id}_recall.table
        """
}

process applyBQSR {

    tag "$sample_id"

    input:
        tuple val(sample_id), path(bam_file)
        tuple val(sample_id), path(recal_table)

    output:
        tuple val(sample_id), path("${sample_id}_recall.bam")

    script:
        """
        gatk ApplyBQSR \
            -R ${params.reference} \
            -I ${bam_file} \
            --bqsr-recal-file ${recal_table} \
            -O ${sample_id}_recall.bam
        """
}
