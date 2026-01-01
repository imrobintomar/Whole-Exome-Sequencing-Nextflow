process fastpQC {
    input:
        tuple val(sample_id), path(reads)
    output:
        tuple val(sample_id), path("*_1_filtered.fastq.gz"), path("*_2_filtered.fastq.gz")
    script:
        """
        fastp -i ${reads[0]} -I ${reads[1]} \\
            --disable_length_filtering \\
            --qualified_quality_phred 30 \\
            -o ${sample_id}_1_filtered.fastq.gz \\
            -O ${sample_id}_2_filtered.fastq.gz \\
            --html ${sample_id}.html \\
            -w 16
        """
}