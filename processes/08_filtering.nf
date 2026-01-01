process extractFilterFields {

    tag "$sample_id"

    input:
        tuple val(sample_id), path(vcf)

    output:
        path("${sample_id}_final_annotated.tsv")

    script:
        """
        # Extract fields and filter variants
        java -jar ${params.snpsift_jar} extractFields -s "," -e "." ${vcf} \\
            CHROM POS ID REF ALT QUAL FILTER DP AC AF AN \\
            ExonicFunc.refGeneWithVer Gene.refGeneWithVer \\
            SIFT_pred Polyphen2_HDIV_pred LRT_pred MutationTaster_pred \\
            FATHMM_pred PROVEAN_pred MetaSVM_pred REVEL_score CADD_phred > \\
            ${sample_id}_annotated.tsv

        # Filter: AF <= max_af (col 10), DP >= min_depth (col 8)
        # Skip header line, then filter based on AF and DP columns
        awk 'NR==1 || (\$10 != "." && \$10 <= ${params.max_af} && \$8 != "." && \$8 >= ${params.min_depth})' \\
            ${sample_id}_annotated.tsv > ${sample_id}_final_annotated.tsv
        """
}