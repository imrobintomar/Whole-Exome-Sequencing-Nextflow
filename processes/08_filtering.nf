process extractFilterFields {

    tag "$sample_id"

    input:
        tuple val(sample_id), path(vcf)

    output:
        path("${sample_id}_final_annotated.tsv")

    script:
        """
        # Extract all ANNOVAR annotated fields without filtering
        java -jar ${params.snpsift_jar} extractFields -s "\t" -e "." ${vcf} \\
            CHROM POS REF ALT QUAL DP \\
            "GEN[0].GT" "GEN[0].AD[1]" "GEN[0].DP" \\
            ExonicFunc.refGeneWithVer AAChange.refGeneWithVer Gene.refGeneWithVer Func.refGeneWithVer \\
            avsnp150 gnomad40_exome_AF ExAC_ALL cosmic84_coding \\
            ALLELEID CLNDISDB CLNDN CLNHGVS CLNREVSTAT CLNSIG CLNVC CLNVCSO \\
            GENEINFO MC ORIGIN RS AF_EXAC AF_ESP CLNSIGCONF AF_TGP CLNVI \\
            CLNDISDBINCL CLNDNINCL CLNSIGINCL DBVARID SCIDISDB SCIDN \\
            ONCDISDB ONCDN \\
            SIFT_pred SIFT4G_pred CADD_phred DANN_rankscore LINSIGHT_rankscore  \\
            Polyphen2_HVAR_pred Polyphen2_HDIV_pred LRT_pred MutationTaster_pred \\
            MutationAssessor_pred FATHMM_pred PROVEAN_pred MetaSVM_pred MetaLR_pred MetaRNN_pred \\
            REVEL_rankscore MutPred_rankscore MVP_rankscore MPC_rankscore \\
            PrimateAI_pred DEOGEN2_pred BayesDel_addAF_pred BayesDel_noAF_pred ClinPred_pred \\
            Interpro_domain GTEx_V8_gene GTEx_V8_tissue \\
            ExcessHet FS MQ QD SOR AC   AF AN   \\
            Ensembl_gene.refGeneWithVer chr.refGeneWithVer Gene_old_names.refGeneWithVer \\
            Gene_other_names.refGeneWithVer Uniprot_acc_HGNC_Uniprot.refGeneWithVer \\
            Uniprot_id_HGNC_Uniprot.refGeneWithVer Entrez_gene_id.refGeneWithVer \\
            CCDS_id.refGeneWithVer Refseq_id.refGeneWithVer ucsc_id.refGeneWithVer \\
            MIM_id.refGeneWithVer OMIM_id.refGeneWithVer Gene_full_name.refGeneWithVer \\
            Pathway_Uniprot.refGeneWithVer Pathway_BioCarta_short.refGeneWithVer \\
            Pathway_BioCarta_full.refGeneWithVer Pathway_ConsensusPathDB.refGeneWithVer \\
            Pathway_KEGG_id.refGeneWithVer Pathway_KEGG_full.refGeneWithVer \\
            Function_description.refGeneWithVer Disease_description.refGeneWithVer \\
            MIM_phenotype_id.refGeneWithVer MIM_disease.refGeneWithVer \\
            Orphanet_disorder_id.refGeneWithVer Orphanet_disorder.refGeneWithVer \\
            Orphanet_association_type.refGeneWithVer Trait_association_GWAS.refGeneWithVer \\
            HPO_id.refGeneWithVer HPO_name.refGeneWithVer \\
            GO_biological_process.refGeneWithVer GO_cellular_component.refGeneWithVer \\
            GO_molecular_function.refGeneWithVer Tissue_specificity_Uniprot.refGeneWithVer \\
            Expression_egenetics.refGeneWithVer Expression_GNF_Atlas.refGeneWithVer \\
            Interactions_IntAct.refGeneWithVer Interactions_BioGRID.refGeneWithVer \\
            Interactions_ConsensusPathDB.refGeneWithVer \\
            P_HI.refGeneWithVer HIPred_score.refGeneWithVer HIPred.refGeneWithVer \\
            GHIS.refGeneWithVer P_rec.refGeneWithVer Known_rec_info.refGeneWithVer \\
            RVIS_EVS.refGeneWithVer RVIS_percentile_EVS.refGeneWithVer \\
            LoF_FDR_ExAC.refGeneWithVer RVIS_ExAC.refGeneWithVer RVIS_percentile_ExAC.refGeneWithVer \\
            ExAC_pLI.refGeneWithVer ExAC_pRec.refGeneWithVer ExAC_pNull.refGeneWithVer \\
            ExAC_nonTCGA_pLI.refGeneWithVer ExAC_nonTCGA_pRec.refGeneWithVer ExAC_nonTCGA_pNull.refGeneWithVer \\
            ExAC_nonpsych_pLI.refGeneWithVer ExAC_nonpsych_pRec.refGeneWithVer ExAC_nonpsych_pNull.refGeneWithVer \\
            gnomAD_pLI.refGeneWithVer gnomAD_pRec.refGeneWithVer gnomAD_pNull.refGeneWithVer \\
            ExAC_del_score.refGeneWithVer ExAC_dup_score.refGeneWithVer ExAC_cnv_score.refGeneWithVer \\
            ExAC_cnv_flag.refGeneWithVer GDI.refGeneWithVer GDI_Phred.refGeneWithVer \\
            Gene_damage_prediction_all_disease_causing_genes.refGeneWithVer \\
            Gene_damage_prediction_all_Mendelian_disease_causing_genes.refGeneWithVer \\
            Gene_damage_prediction_Mendelian_AD_disease_causing_genes.refGeneWithVer \\
            Gene_damage_prediction_Mendelian_AR_disease_causing_genes.refGeneWithVer \\
            Gene_damage_prediction_all_PID_disease_causing_genes.refGeneWithVer \\
            Gene_damage_prediction_PID_AD_disease_causing_genes.refGeneWithVer \\
            Gene_damage_prediction_PID_AR_disease_causing_genes.refGeneWithVer \\
            Gene_damage_prediction_all_cancer_disease_causing_genes.refGeneWithVer \\
            Gene_damage_prediction_cancer_recessive_disease_causing_genes.refGeneWithVer \\
            Gene_damage_prediction_cancer_dominant_disease_causing_genes.refGeneWithVer \\
            LoFtool_score.refGeneWithVer \\
            SORVA_LOF_MAF0.005_HetOrHom.refGeneWithVer SORVA_LOF_MAF0.005_HomOrCompoundHet.refGeneWithVer \\
            SORVA_LOF_MAF0.001_HetOrHom.refGeneWithVer SORVA_LOF_MAF0.001_HomOrCompoundHet.refGeneWithVer \\
            SORVA_LOForMissense_MAF0.005_HetOrHom.refGeneWithVer \\
            SORVA_LOForMissense_MAF0.005_HomOrCompoundHet.refGeneWithVer \\
            SORVA_LOForMissense_MAF0.001_HetOrHom.refGeneWithVer \\
            SORVA_LOForMissense_MAF0.001_HomOrCompoundHet.refGeneWithVer \\
            ZFIN_zebrafish_gene.refGeneWithVer ZFIN_zebrafish_structure.refGeneWithVer \\
            FILTER "ANN[0].ALLELE" "ANN[0].EFFECT" "ANN[0].IMPACT" "ANN[0].GENE" "ANN[0].RANK" \\
            "EFF[0].CODON" "ANN[0].FEATUREID" "ANN[0].GENEID" "ANN[0].HGVS_C" "ANN[0].HGVS_P" "ANN[0].BIOTYPE" \\
            gnomad40_exome_AF_raw gnomad40_exome_AF_XX gnomad40_exome_AF_XY gnomad40_exome_AF_grpmax \\
            gnomad40_exome_faf95 gnomad40_exome_faf99 gnomad40_exome_fafmax_faf99_max \\
            gnomad40_exome_AF_afr gnomad40_exome_AF_amr gnomad40_exome_AF_asj gnomad40_exome_AF_eas \\
            gnomad40_exome_AF_fin gnomad40_exome_AF_mid gnomad40_exome_AF_nfe \\
            gnomad40_exome_AF_remaining gnomad40_exome_AF_sas \\
            ExAC_AFR ExAC_AMR ExAC_EAS ExAC_FIN ExAC_NFE ExAC_OTH ExAC_SAS \\
            SIFT_score SIFT_converted_rankscore SIFT4G_score SIFT4G_converted_rankscore \\
            Polyphen2_HDIV_score Polyphen2_HDIV_rankscore Polyphen2_HVAR_score Polyphen2_HVAR_rankscore \\
            LRT_score LRT_converted_rankscore MutationTaster_score MutationTaster_converted_rankscore \\
            MutationAssessor_score MutationAssessor_rankscore FATHMM_score FATHMM_converted_rankscore \\
            PROVEAN_score PROVEAN_converted_rankscore VEST4_score VEST4_rankscore \\
            MetaSVM_score MetaSVM_rankscore MetaLR_score MetaLR_rankscore \\
            MetaRNN_score MetaRNN_rankscore REVEL_score MutPred_score MVP_score MPC_score \\
            PrimateAI_score PrimateAI_rankscore DEOGEN2_score DEOGEN2_rankscore \\
            BayesDel_addAF_score BayesDel_addAF_rankscore BayesDel_noAF_score BayesDel_noAF_rankscore \\
            ClinPred_score ClinPred_rankscore Aloft_pred Aloft_Confidence \\
            CADD_raw CADD_raw_rankscore \\
            > ${sample_id}_extracted.tsv

        # Create UniqueID column (chr:pos:ref:alt) and insert after POS column
        awk 'BEGIN {FS=OFS="\\t"}
        NR==1 {
            # Header line: insert "UniqueID" after POS (column 2)
            print \$1, \$2, "UniqueID", \$3, \$4, \$5, \$6, \$7, \$8, \$9;
            for(i=10; i<=NF; i++) printf "%s%s", OFS, \$i;
            print "";
            next;
        }
        {
            # Data lines: create UniqueID as chr:pos:ref:alt
            uniqueID = \$1 ":" \$2 ":" \$3 ":" \$4;
            print \$1, \$2, uniqueID, \$3, \$4, \$5, \$6, \$7, \$8, \$9;
            for(i=10; i<=NF; i++) printf "%s%s", OFS, \$i;
            print "";
        }' ${sample_id}_extracted.tsv > ${sample_id}_final_annotated.tsv
        """
}