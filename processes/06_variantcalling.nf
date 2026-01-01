process haplotypeCaller {

    tag "$sample_id"

    input:
        tuple val(sample_id), path(bam_file), path(bam_index)

    output:
        tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.tbi")

    script:
        """
        gatk HaplotypeCaller -R ${params.reference} \\
            -I ${bam_file} \\
            -O ${sample_id}.vcf.gz \\
            --native-pair-hmm-threads ${task.cpus}
        """
}