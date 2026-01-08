process annovarAnnotate {

    tag "$sample_id"

    input:
        tuple val(sample_id), path(vcf), path(vcf_index)

    output:
        tuple val(sample_id), path("${sample_id}.annovar.hg38_multianno.txt")

    script:
        """
        # ANNOVAR as single source of truth - outputs TXT format directly
        ${params.annovar_dir}/table_annovar.pl ${vcf} ${params.annovar_db} \\
            --buildver hg38 \\
            --out ${sample_id}.annovar \\
            --remove \\
            --protocol refGeneWithVer,dbnsfp42a,clinvar_20240416,gnomad40_exome,avsnp150,cosmic84_coding,exac03 \\
            --operation gx,f,f,f,f,f,f \\
            -nastring . \\
            --polish \\
            --xreffile ${params.annovar_xreffile} \\
            --thread ${task.cpus}
        """
}