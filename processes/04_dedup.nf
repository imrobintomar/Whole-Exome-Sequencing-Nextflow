/*
 * Duplicate marking and post-deduplication sorting
 */

process markDuplicates {

    tag "$sample_id"

    input:
        tuple val(sample_id), path(bam_file), path(bam_index)

    output:
        tuple val(sample_id), path("${sample_id}_markdup.bam"), path("${sample_id}_markdup.metrics.txt")

    script:
        // Calculate safe Java heap size (90% of allocated memory)
        def avail_mem_mb = (task.memory.toMega() * 0.9).toInteger()
        """
        mkdir -p tmp_${sample_id}

        gatk --java-options "-Xmx${avail_mem_mb}m -XX:+UseParallelGC -Djava.io.tmpdir=./tmp_${sample_id}" \
            MarkDuplicates \
            -I ${bam_file} \
            -O ${sample_id}_markdup.bam \
            -M ${sample_id}_markdup.metrics.txt \
            --TMP_DIR ./tmp_${sample_id}

        # Cleanup
        rm -rf tmp_${sample_id}
        """
}

process sortSamPostDedup {

    tag "$sample_id"

    input:
        tuple val(sample_id), path(bam_file), path(metrics)

    output:
        tuple val(sample_id), path("${sample_id}_markdup.sorted.bam"), path("${sample_id}_markdup.sorted.bam.bai")

    script:
        // Calculate safe Java heap size (90% of allocated memory)
        def avail_mem_mb = (task.memory.toMega() * 0.9).toInteger()
        """
        mkdir -p tmp_${sample_id}

        gatk --java-options "-Xmx${avail_mem_mb}m -XX:+UseParallelGC -Djava.io.tmpdir=./tmp_${sample_id}" \
            SortSam \
            -I ${bam_file} \
            -O ${sample_id}_markdup.sorted.bam \
            -SO coordinate \
            --CREATE_INDEX false \
            --TMP_DIR ./tmp_${sample_id}

        # Use samtools for reliable index naming
        samtools index ${sample_id}_markdup.sorted.bam

        # Cleanup
        rm -rf tmp_${sample_id}
        """
}
