process bwaMem {

    tag "${sample_id}"

    input:
        tuple val(sample_id), path(read1), path(read2)

    output:
        tuple val(sample_id), path("${sample_id}.bam")

    script:
        """
        # BWA-MEM alignment with proper flags:
        # -M: Mark shorter split hits as secondary (required for GATK MarkDuplicates)
        # -Y: Use soft clipping for supplementary alignments
        # -K: Process 100M input bases per batch for deterministic output
        bwa mem -M -Y -K 100000000 -t ${task.cpus} \\
            -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tLB:${sample_id}\\tPL:ILLUMINA" \\
            ${params.reference} ${read1} ${read2} | \\
            samtools view -@ 4 -bh -o ${sample_id}.bam -
        """
}
