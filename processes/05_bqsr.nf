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
        // Calculate safe Java heap size (90% of allocated memory)
        def avail_mem_mb = (task.memory.toMega() * 0.9).toInteger()
        // Add intervals if provided (critical for exome sequencing)
        def intervals_arg = params.intervals ? "-L ${params.intervals}" : ""

        """
        mkdir -p tmp_${sample_id}

        gatk --java-options "-Xmx${avail_mem_mb}m -XX:+UseParallelGC -Djava.io.tmpdir=./tmp_${sample_id}" \
            BaseRecalibrator \
            -R ${params.reference} \
            -I ${bam_file} \
            ${intervals_arg} \
            ${known_sites_str} \
            -O ${sample_id}_recall.table \
            --tmp-dir ./tmp_${sample_id}

        # Cleanup
        rm -rf tmp_${sample_id}
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
        // Calculate safe Java heap size (90% of allocated memory)
        def avail_mem_mb = (task.memory.toMega() * 0.9).toInteger()
        """
        mkdir -p tmp_${sample_id}

        gatk --java-options "-Xmx${avail_mem_mb}m -XX:+UseParallelGC -Djava.io.tmpdir=./tmp_${sample_id}" \
            ApplyBQSR \
            -R ${params.reference} \
            -I ${bam_file} \
            --bqsr-recal-file ${recal_table} \
            -O ${sample_id}_recall.bam \
            --create-output-bam-index false \
            --tmp-dir ./tmp_${sample_id}

        # Use samtools for reliable index naming
        samtools index ${sample_id}_recall.bam

        # Cleanup
        rm -rf tmp_${sample_id}
        """
}
