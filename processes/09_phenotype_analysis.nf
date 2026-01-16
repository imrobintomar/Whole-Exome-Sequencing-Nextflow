process exomiserPhenotypeAnalysis {

    tag "$sample_id"
    publishDir "${params.output_dir}/Germline_VCF/${sample_id}", mode: 'copy'

    input:
        tuple val(sample_id), path(vcf), path(annovar_txt), val(hpo_terms)

    output:
        tuple val(sample_id), path("${sample_id}.hg38_multianno.hpo.txt"), emit: augmented_annovar
        tuple val(sample_id), path("${sample_id}_exomiser.tsv"), emit: exomiser_output
        tuple val(sample_id), path("${sample_id}_hpo_gene_rank.png"), emit: gene_rank_plot
        tuple val(sample_id), path("${sample_id}_hpo_score_dist.png"), emit: score_dist_plot

    script:
        def hpo_list = hpo_terms instanceof List ? hpo_terms.join(',') : hpo_terms
        """
        #!/usr/bin/env python3
        import os
        import sys
        import yaml
        import pandas as pd
        import subprocess
        import plotly.graph_objects as go
        import plotly.express as px
        from pathlib import Path

        # Configuration
        sample_id = "${sample_id}"
        vcf_path = "${vcf}"
        annovar_path = "${annovar_txt}"
        hpo_terms = "${hpo_list}".split(',')
        exomiser_dir = "${params.exomiser_dir}"
        exomiser_data = "${params.exomiser_data_dir}"

        # Generate Exomiser YAML
        exomiser_config = {
            'analysis': {
                'vcf': str(Path(vcf_path).absolute()),
                'genomeAssembly': 'hg38',
                'hpoIds': hpo_terms,
                'inheritanceModes': {
                    'AUTOSOMAL_DOMINANT': 1.0,
                    'AUTOSOMAL_RECESSIVE': 1.0,
                    'X_RECESSIVE': 1.0,
                    'X_DOMINANT': 1.0,
                    'MITOCHONDRIAL': 1.0
                },
                'analysisMode': 'PASS_ONLY',
                'frequencySources': ['THOUSAND_GENOMES', 'TOPMED', 'GNOMAD_E_AFR', 'GNOMAD_E_AMR', 'GNOMAD_E_EAS', 'GNOMAD_E_FIN', 'GNOMAD_E_NFE', 'GNOMAD_E_SAS'],
                'pathogenicitySources': ['REVEL', 'MVP'],
                'steps': [
                    {'variantEffectFilter': {'remove': ['SYNONYMOUS_VARIANT', 'UPSTREAM_GENE_VARIANT', 'DOWNSTREAM_GENE_VARIANT', 'INTERGENIC_VARIANT', 'INTRON_VARIANT']}},
                    {'frequencyFilter': {'maxFrequency': 0.01}},
                    {'pathogenicityFilter': {'keepNonPathogenic': True}},
                    {'inheritanceFilter': {}},
                    {'omimPrioritiser': {}},
                    {'hiPhivePrioritiser': {}}
                ]
            },
            'outputOptions': {
                'outputFormats': ['TSV'],
                'outputPrefix': sample_id
            }
        }

        yaml_path = f"{sample_id}_exomiser.yml"
        with open(yaml_path, 'w') as f:
            yaml.dump(exomiser_config, f)

        # Run Exomiser
        print(f"Running Exomiser with HPO terms: {hpo_terms}")

        # Find the correct JAR file (handles version-suffixed names)
        import glob
        jar_files = glob.glob(f'{exomiser_dir}/exomiser-cli*.jar')
        if not jar_files:
            raise FileNotFoundError(f"No Exomiser JAR found in {exomiser_dir}")
        exomiser_jar = jar_files[0]
        print(f"Using Exomiser JAR: {exomiser_jar}")

        exomiser_cmd = [
            'java', '-Xms2g', '-Xmx4g',
            '-jar', exomiser_jar,
            '--analysis', yaml_path,
            '--assembly', 'hg38'
        ]

        result = subprocess.run(exomiser_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"WARNING: Exomiser failed: {result.stderr}")
            # Create empty output on failure
            exomiser_tsv = f"{sample_id}_exomiser.tsv"
            with open(exomiser_tsv, 'w') as f:
                f.write("#Exomiser_failed\\n")
        else:
            # Find Exomiser TSV output
            exomiser_tsv = f"{sample_id}.variants.tsv"
            if not os.path.exists(exomiser_tsv):
                print(f"WARNING: Exomiser output not found at {exomiser_tsv}")
                with open(f"{sample_id}_exomiser.tsv", 'w') as f:
                    f.write("#Exomiser_output_missing\\n")
                exomiser_tsv = f"{sample_id}_exomiser.tsv"

        # Parse Exomiser results
        def parse_exomiser(tsv_path):
            if not os.path.exists(tsv_path) or os.path.getsize(tsv_path) == 0:
                return pd.DataFrame()

            try:
                df = pd.read_csv(tsv_path, sep='\\t', comment='#')
                if df.empty:
                    return pd.DataFrame()

                # Extract key columns
                exomiser_data = {}
                for _, row in df.iterrows():
                    chrom = str(row.get('CHROM', row.get('#CHROM', '')))
                    pos = str(row.get('POS', ''))
                    ref = str(row.get('REF', ''))
                    alt = str(row.get('ALT', ''))
                    key = f"{chrom}:{pos}:{ref}:{alt}"

                    exomiser_data[key] = {
                        'HPO_GENE': row.get('GENE_SYMBOL', '.'),
                        'HPO_DISEASE': row.get('CONTRIBUTING_DISEASE', '.'),
                        'HPO_SCORE': row.get('EXOMISER_GENE_COMBINED_SCORE', 0.0),
                        'HPO_PHENO_SCORE': row.get('EXOMISER_GENE_PHENO_SCORE', 0.0),
                        'HPO_VARIANT_SCORE': row.get('EXOMISER_VARIANT_SCORE', 0.0),
                        'HPO_INHERITANCE': row.get('CONTRIBUTING_MOI', '.')
                    }
                return exomiser_data
            except Exception as e:
                print(f"Error parsing Exomiser output: {e}")
                return {}

        exomiser_results = parse_exomiser(exomiser_tsv)

        # Load ANNOVAR file
        print(f"Loading ANNOVAR file: {annovar_path}")
        annovar_df = pd.read_csv(annovar_path, sep='\\t', encoding='utf-8', on_bad_lines='skip')

        # Merge Exomiser results
        new_cols = ['HPO_GENE', 'HPO_DISEASE', 'HPO_SCORE', 'HPO_PHENO_SCORE', 'HPO_VARIANT_SCORE', 'HPO_INHERITANCE']
        for col in new_cols:
            annovar_df[col] = '.'

        if exomiser_results:
            print(f"Merging {len(exomiser_results)} Exomiser variants with ANNOVAR")
            for idx, row in annovar_df.iterrows():
                chrom = str(row.get('Chr', row.get('CHROM', '')))
                pos = str(row.get('Start', row.get('POS', '')))
                ref = str(row.get('Ref', row.get('REF', '')))
                alt = str(row.get('Alt', row.get('ALT', '')))
                key = f"{chrom}:{pos}:{ref}:{alt}"

                if key in exomiser_results:
                    for col, val in exomiser_results[key].items():
                        annovar_df.at[idx, col] = val
        else:
            print("No Exomiser results to merge")

        # Save augmented ANNOVAR file
        output_path = f"{sample_id}.hg38_multianno.hpo.txt"
        annovar_df.to_csv(output_path, sep='\\t', index=False)
        print(f"Saved augmented ANNOVAR: {output_path}")

        # Generate visualizations
        def safe_float(val):
            try:
                return float(val) if pd.notna(val) and val != '.' else 0.0
            except:
                return 0.0

        # Gene ranking plot
        gene_scores = {}
        for idx, row in annovar_df.iterrows():
            gene = row['HPO_GENE']
            score = safe_float(row['HPO_SCORE'])
            if gene != '.' and score > 0:
                if gene not in gene_scores:
                    gene_scores[gene] = score
                else:
                    gene_scores[gene] = max(gene_scores[gene], score)

        if gene_scores:
            top_genes = sorted(gene_scores.items(), key=lambda x: x[1], reverse=True)[:20]
            genes, scores = zip(*top_genes)

            fig = go.Figure(data=[go.Bar(x=list(scores), y=list(genes), orientation='h', marker_color='#06b6d4')])
            fig.update_layout(
                title=f'Top 20 Genes by HPO Score - {sample_id}',
                xaxis_title='Exomiser Combined Score',
                yaxis_title='Gene',
                height=600,
                template='plotly_white'
            )
            fig.write_image(f"{sample_id}_hpo_gene_rank.png", width=1200, height=800)
        else:
            # Empty plot
            fig = go.Figure()
            fig.add_annotation(text="No phenotype-matched genes", xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False)
            fig.write_image(f"{sample_id}_hpo_gene_rank.png", width=800, height=600)

        # Score distribution plot
        scores = [safe_float(row['HPO_SCORE']) for idx, row in annovar_df.iterrows() if safe_float(row['HPO_SCORE']) > 0]

        if scores:
            fig = go.Figure(data=[go.Histogram(x=scores, nbinsx=30, marker_color='#14b8a6')])
            fig.update_layout(
                title=f'HPO Score Distribution - {sample_id}',
                xaxis_title='Exomiser Combined Score',
                yaxis_title='Variant Count',
                template='plotly_white'
            )
            fig.write_image(f"{sample_id}_hpo_score_dist.png", width=1000, height=600)
        else:
            fig = go.Figure()
            fig.add_annotation(text="No variants with phenotype scores", xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False)
            fig.write_image(f"{sample_id}_hpo_score_dist.png", width=800, height=600)

        print("Phenotype analysis complete")
        """
}
