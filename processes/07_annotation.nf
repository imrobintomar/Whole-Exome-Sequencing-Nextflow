process annovarAnnotate {

    tag "$sample_id"

    input:
        tuple val(sample_id), path(vcf), path(vcf_index)

    output:
        tuple val(sample_id), path("${sample_id}.annovar.hg38_multianno.txt")

    script:
        """
        # Uncompress VCF if it's gzipped (ANNOVAR doesn't handle .vcf.gz properly)
        if [[ ${vcf} == *.gz ]]; then
            echo "Uncompressing gzipped VCF: ${vcf}"
            gunzip -c ${vcf} > ${sample_id}_input.vcf
            INPUT_VCF="${sample_id}_input.vcf"
        else
            INPUT_VCF="${vcf}"
        fi

        # ANNOVAR as single source of truth - outputs TXT format directly
        ${params.annovar_dir}/table_annovar.pl \$INPUT_VCF ${params.annovar_db} \\
            --buildver hg38 \\
            --out ${sample_id}.annovar \\
            --remove \\
            --protocol refGeneWithVer,dbnsfp42a,clinvar_20240416,gnomad40_exome,avsnp150,cosmic84_coding,exac03 \\
            --operation gx,f,f,f,f,f,f \\
            -nastring . \\
            --polish \\
            --xreffile ${params.annovar_xreffile} \\
            --thread ${task.cpus}

        # Validate output is not corrupted
        OUTPUT_FILE="${sample_id}.annovar.hg38_multianno.txt"
        if [ ! -f "\$OUTPUT_FILE" ]; then
            echo "ERROR: ANNOVAR output file not created!"
            exit 1
        fi

        # Check if output is valid text (not binary)
        if ! file "\$OUTPUT_FILE" | grep -qi "text"; then
            echo "ERROR: ANNOVAR output is corrupted (binary data detected)!"
            echo "File type: \$(file \$OUTPUT_FILE)"
            exit 1
        fi

        # Check for excessive invalid input
        INVALID_FILE="${sample_id}.annovar.invalid_input"
        if [ -f "\$INVALID_FILE" ]; then
            INVALID_COUNT=\$(wc -l < "\$INVALID_FILE")
            TOTAL_VARIANTS=\$(grep -v "^#" \$INPUT_VCF | wc -l)
            if [ \$INVALID_COUNT -eq \$TOTAL_VARIANTS ] && [ \$TOTAL_VARIANTS -gt 0 ]; then
                echo "ERROR: ALL variants marked as invalid by ANNOVAR!"
                echo "This indicates a format or compression issue."
                exit 1
            fi
            if [ \$INVALID_COUNT -gt 0 ]; then
                echo "WARNING: \$INVALID_COUNT variants had invalid input format"
            fi
        fi

        # Clean up temporary uncompressed VCF
        if [[ ${vcf} == *.gz ]]; then
            rm -f ${sample_id}_input.vcf
        fi

        echo "ANNOVAR annotation completed successfully"
        """
}