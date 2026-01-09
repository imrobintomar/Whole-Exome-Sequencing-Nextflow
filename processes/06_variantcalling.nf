process haplotypeCaller {

    tag "$sample_id"

    input:
        tuple val(sample_id), path(bam_file), path(bam_index)

    output:
        tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.tbi")

    script:
        // Calculate safe Java heap size (90% of allocated memory)
        def avail_mem_mb = (task.memory.toMega() * 0.9).toInteger()
        // Add intervals if provided (CRITICAL for exome sequencing - 60x speedup)
        def intervals_arg = params.intervals ? "-L ${params.intervals}" : ""
        """
        mkdir -p tmp_${sample_id}

        gatk --java-options "-Xmx${avail_mem_mb}m -XX:+UseParallelGC -Djava.io.tmpdir=./tmp_${sample_id}" \\
            HaplotypeCaller \\
            -R ${params.reference} \\
            -I ${bam_file} \\
            ${intervals_arg} \\
            -O ${sample_id}.vcf.gz \\
            --native-pair-hmm-threads ${task.cpus} \\
            --tmp-dir ./tmp_${sample_id}

        # Cleanup
        rm -rf tmp_${sample_id}
        """
}