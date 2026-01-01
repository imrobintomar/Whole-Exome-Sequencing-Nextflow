/*
 * Duplicate marking and post-deduplication sorting
 */

process markDuplicates {

    tag "$sample_id"

    input:
        tuple val(sample_id), path(bam_file)

    output:
        tuple val(sample_id), path("${sample_id}_markdup.bam")

    script:
        """
        gatk MarkDuplicates \
            -I ${bam_file} \
            -O ${sample_id}_markdup.bam \
            -M ${sample_id}_markdup.metrics.txt
        """
}

process sortSamPostDedup {

    tag "$sample_id"

    input:
        tuple val(sample_id), path(bam_file)

    output:
        tuple val(sample_id), path("${sample_id}_markdup.sorted.bam"), path("${sample_id}_markdup.sorted.bam.bai")

    script:
        """
        gatk SortSam \
            -I ${bam_file} \
            -O ${sample_id}_markdup.sorted.bam \
            -SO coordinate \
            --CREATE_INDEX true
        """
}
