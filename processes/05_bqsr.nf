process baseRecalibrator {

    tag "$sample_id"

    input:
        tuple val(sample_id), path(bam_file), path(bam_index)

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
        tuple val(sample_id), path(bam_file), path(bam_index)
        tuple val(sample_id), path(recal_table)

    output:
        tuple val(sample_id), path("${sample_id}_recall.bam"), path("${sample_id}_recall.bam.bai")

    script:
        """
        gatk ApplyBQSR \
            -R ${params.reference} \
            -I ${bam_file} \
            --bqsr-recal-file ${recal_table} \
            -O ${sample_id}_recall.bam \
            --create-output-bam-index true

        # Ensure index has correct name
        if [ -f ${sample_id}_recall.bai ]; then
            mv ${sample_id}_recall.bai ${sample_id}_recall.bam.bai
        fi
        """
}
