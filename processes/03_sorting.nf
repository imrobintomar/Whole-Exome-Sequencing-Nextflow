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
        // Calculate safe Java heap size (90% of allocated memory)
        def avail_mem_mb = (task.memory.toMega() * 0.9).toInteger()
        """
        mkdir -p tmp_${sample_id}

        gatk --java-options "-Xmx${avail_mem_mb}m -XX:+UseParallelGC -Djava.io.tmpdir=./tmp_${sample_id}" \
            SortSam \
            -I ${bam_in} \
            -O ${sample_id}.sorted.bam \
            -SO coordinate \
            --CREATE_INDEX false \
            --TMP_DIR ./tmp_${sample_id}

        # Use samtools for reliable index naming
        samtools index ${sample_id}.sorted.bam

        # Cleanup
        rm -rf tmp_${sample_id}
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
