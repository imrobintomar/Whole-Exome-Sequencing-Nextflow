process snpsiftAnnotate1000G {

    tag "$sample_id"

    input:
        tuple val(sample_id), path(vcf), path(vcf_index)

    output:
        tuple val(sample_id), path("${sample_id}_1000genome.vcf.gz")

    script:
        """
        # Decompress input VCF if needed
        if [[ ${vcf} == *.gz ]]; then
            gunzip -c ${vcf} > temp_input.vcf
            vcf_temp="temp_input.vcf"
        else
            vcf_temp="${vcf}"
        fi

        # Annotate with 1000 Genomes data
        for chr in {1..22} X; do
            vcf_file="${params.thousand_genomes_dir}/filtered_chr\${chr}.vcf"
            if [[ -f "\${vcf_file}" ]] || [[ -f "\${vcf_file}.gz" ]]; then
                java -jar ${params.snpsift_jar} annotate \\
                    \${vcf_file} \\
                    \${vcf_temp} > \${vcf_temp}.tmp && \\
                mv \${vcf_temp}.tmp \${vcf_temp}
            else
                echo "Warning: 1000G file for chr\${chr} not found, skipping..."
            fi
        done

        # Compress and index output
        bgzip -c \${vcf_temp} > ${sample_id}_1000genome.vcf.gz
        tabix -p vcf ${sample_id}_1000genome.vcf.gz
        """
}

process annovarAnnotate {

    tag "$sample_id"

    input:
        tuple val(sample_id), path(vcf)

    output:
        tuple val(sample_id), path("${sample_id}.annovar.hg38_multianno.vcf")

    script:
        """
        ${params.annovar_dir}/table_annovar.pl ${vcf} ${params.annovar_db} \\
            --buildver hg38 \\
            --out ${sample_id}.annovar \\
            --remove \\
            --protocol refGeneWithVer,dbnsfp42a,clinvar_20240416,gnomad40_exome,avsnp150,cosmic84_coding,exac03 \\
            --operation gx,f,f,f,f,f,f \\
            -nastring . \\
            --polish \\
            --xreffile ${params.annovar_xreffile} \\
            --vcfinput \\
            --thread ${task.cpus}
        """
}