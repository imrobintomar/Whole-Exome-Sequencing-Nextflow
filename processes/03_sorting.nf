/*
 * Sorting & alignment QC
 * Used ONLY after alignment
 */

process sortSam {

    tag "$sample_id"

    input:
        tuple val(sample_id), path(bam_in)

    output:
        tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai")

    script:
        """
        gatk SortSam \
            -I ${bam_in} \
            -O ${sample_id}.sorted.bam \
            -SO coordinate \
            --CREATE_INDEX true

        mv ${sample_id}.sorted.bai ${sample_id}.sorted.bam.bai
        """
}

process flagstat {

    tag "$sample_id"

    input:
        tuple val(sample_id), path(bam_file), path(bam_index)

    output:
        path("${sample_id}.Stat.txt")

    script:
        """
        samtools flagstat ${bam_file} > ${sample_id}.Stat.txt
        """
}
