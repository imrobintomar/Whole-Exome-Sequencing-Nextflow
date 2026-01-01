process bwaMem {

    tag "${sample_id}"

    input:
        tuple val(sample_id), path(read1), path(read2)

    output:
        tuple val(sample_id), path("${sample_id}.bam")

    script:
        """
        bwa mem -Y -K 100000000 -t ${task.cpus} \\
            -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tLB:${sample_id}\\tPL:ILLUMINA" \\
            ${params.reference} ${read1} ${read2} | \\
            samtools view -bh -o ${sample_id}.bam -
        """
}
