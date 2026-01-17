"""
Variant Analysis Module
Analyzes filtered TSV files to extract clinical/QC metrics
"""

import pandas as pd
from typing import Dict, List, Any
from pathlib import Path


class VariantAnalyzer:
    """Analyzes variant TSV files and extracts metrics"""

    # Chromosome order for plotting
    CHROMOSOMES = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY', 'chrM']

    # Chromosome lengths (hg38) in Mb
    CHR_LENGTHS = {
        'chr1': 248.96, 'chr2': 242.19, 'chr3': 198.30, 'chr4': 190.21,
        'chr5': 181.54, 'chr6': 170.81, 'chr7': 159.35, 'chr8': 145.14,
        'chr9': 138.39, 'chr10': 133.80, 'chr11': 135.09, 'chr12': 133.28,
        'chr13': 114.36, 'chr14': 107.04, 'chr15': 101.99, 'chr16': 90.34,
        'chr17': 83.26, 'chr18': 80.37, 'chr19': 58.62, 'chr20': 64.44,
        'chr21': 46.71, 'chr22': 50.82, 'chrX': 156.04, 'chrY': 57.23, 'chrM': 0.016
    }

    def __init__(self, tsv_path: str):
        """Initialize with TSV file path"""
        self.tsv_path = Path(tsv_path)
        self.df = None
        self._load_data()

    def _load_data(self):
        """Load TSV file into pandas DataFrame"""
        if not self.tsv_path.exists():
            raise FileNotFoundError(f"TSV file not found: {self.tsv_path}")

        # Handle encoding and parsing issues from pipeline output
        try:
            self.df = pd.read_csv(self.tsv_path, sep='\t', low_memory=False, encoding='utf-8', on_bad_lines='skip')
        except UnicodeDecodeError:
            self.df = pd.read_csv(self.tsv_path, sep='\t', low_memory=False, encoding='latin1', on_bad_lines='skip')
        except Exception:
            # Last resort: use Python engine with error handling
            self.df = pd.read_csv(self.tsv_path, sep='\t', low_memory=False, encoding='latin1', on_bad_lines='skip', engine='python')

    def _is_snv(self, ref: str, alt: str) -> bool:
        """Check if variant is SNV (both REF and ALT are single nucleotides)"""
        if pd.isna(ref) or pd.isna(alt):
            return False
        ref = str(ref).strip()
        alt = str(alt).strip()
        return len(ref) == 1 and len(alt) == 1

    def _is_indel(self, ref: str, alt: str) -> bool:
        """Check if variant is INDEL (length difference between REF and ALT)"""
        if pd.isna(ref) or pd.isna(alt):
            return False
        ref = str(ref).strip()
        alt = str(alt).strip()
        return len(ref) != len(alt)

    def get_headline_metrics(self) -> Dict[str, int]:
        """
        Calculate headline metrics:
        - Total variants
        - Total unique genes
        - SNVs
        - INDELs
        - High-confidence (PASS) variants
        """
        total_variants = len(self.df)

        # Count unique genes (handle multiple genes per variant)
        gene_col = None
        for col in ['Gene.refGeneWithVer', 'Gene.refGene', 'Gene']:
            if col in self.df.columns:
                gene_col = col
                break

        if gene_col:
            # Split comma-separated genes and count unique
            all_genes = set()
            for genes_str in self.df[gene_col].dropna():
                genes = str(genes_str).split(',')
                all_genes.update(g.strip() for g in genes if g.strip())
            total_genes = len(all_genes)
        else:
            total_genes = 0

        # Count SNVs and INDELs
        ref_col = 'Ref' if 'Ref' in self.df.columns else 'REF'
        alt_col = 'Alt' if 'Alt' in self.df.columns else 'ALT'

        snvs = 0
        indels = 0
        if ref_col in self.df.columns and alt_col in self.df.columns:
            for _, row in self.df.iterrows():
                if self._is_snv(row[ref_col], row[alt_col]):
                    snvs += 1
                elif self._is_indel(row[ref_col], row[alt_col]):
                    indels += 1

        # Count PASS variants (both "PASS" and "." are considered high confidence)
        filter_col = 'FILTER' if 'FILTER' in self.df.columns else 'Filter'
        high_confidence = 0
        if filter_col in self.df.columns:
            high_confidence = len(self.df[self.df[filter_col].isin(['PASS', '.'])])

        return {
            'total_variants': total_variants,
            'total_genes': total_genes,
            'snvs': snvs,
            'indels': indels,
            'high_confidence': high_confidence
        }

    def get_chromosome_distribution(self, pass_only: bool = True, normalize: bool = False) -> List[Dict[str, Any]]:
        """
        Get variant counts per chromosome

        Args:
            pass_only: Only count PASS variants
            normalize: Normalize by chromosome length (variants per Mb)

        Returns:
            List of dicts with chromosome, count, and type (autosome/sex)
        """
        chr_col = None
        for col in ['Chr', 'CHROM', '#CHROM', 'Chromosome']:
            if col in self.df.columns:
                chr_col = col
                break

        if not chr_col:
            return []

        df = self.df.copy()

        # Filter for PASS only (both "PASS" and "." are considered high confidence)
        if pass_only:
            filter_col = 'FILTER' if 'FILTER' in self.df.columns else 'Filter'
            if filter_col in df.columns:
                df = df[df[filter_col].isin(['PASS', '.'])]

        # Count variants per chromosome
        chr_counts = df[chr_col].value_counts().to_dict()

        result = []
        for chr_name in self.CHROMOSOMES:
            count = chr_counts.get(chr_name, 0)

            # Normalize if requested
            if normalize and chr_name in self.CHR_LENGTHS:
                count = count / self.CHR_LENGTHS[chr_name]

            # Determine type
            chr_type = 'sex' if chr_name in ['chrX', 'chrY'] else 'autosome'

            result.append({
                'chromosome': chr_name,
                'count': count,
                'type': chr_type
            })

        return result

    def get_gene_distribution(self, top_n: int = 20, protein_altering_only: bool = False) -> List[Dict[str, Any]]:
        """
        Get top genes by variant count

        Args:
            top_n: Number of top genes to return
            protein_altering_only: Only count protein-altering variants

        Returns:
            List of dicts with gene name and count
        """
        gene_col = None
        for col in ['Gene.refGeneWithVer', 'Gene.refGene', 'Gene']:
            if col in self.df.columns:
                gene_col = col
                break

        if not gene_col:
            return []

        df = self.df.copy()

        # Filter for protein-altering variants
        if protein_altering_only:
            # Try both ExonicFunc column name variants
            func_col = None
            for col in ['ExonicFunc.refGeneWithVer', 'ExonicFunc.refGene', 'ExonicFunc']:
                if col in df.columns:
                    func_col = col
                    break

            if func_col:
                # ANNOVAR ExonicFunc values that indicate protein-altering:
                # - nonsynonymous SNV (amino acid change)
                # - stopgain (premature stop codon)
                # - stoploss (loss of stop codon)
                # - frameshift insertion/deletion
                # - nonframeshift insertion/deletion (in-frame indel)
                # - splicing (affects splice sites)
                protein_altering = [
                    'nonsynonymous',  # matches "nonsynonymous SNV"
                    'stopgain',
                    'stoploss',
                    'frameshift',     # matches both frameshift insertion and deletion
                    'nonframeshift',  # matches nonframeshift insertion/deletion
                    'splicing',
                    'missense',       # alternative term
                    'nonsense'        # alternative term
                ]
                df = df[df[func_col].str.contains('|'.join(protein_altering), case=False, na=False)]

        # Count variants per gene
        gene_counts = {}
        for genes_str in df[gene_col].dropna():
            genes = str(genes_str).split(',')
            for gene in genes:
                gene = gene.strip()
                if gene and gene != '.':
                    gene_counts[gene] = gene_counts.get(gene, 0) + 1

        # Sort and get top N
        sorted_genes = sorted(gene_counts.items(), key=lambda x: x[1], reverse=True)[:top_n]

        return [{'gene': gene, 'count': count} for gene, count in sorted_genes]

    def get_functional_impact(self) -> Dict[str, Any]:
        """
        Get functional impact summary

        Returns:
            Dict with categories and their counts
        """
        # Try to find the functional annotation column
        func_col = None
        for col in ['Func.refGeneWithVer', 'Func.refGene', 'Func']:
            if col in self.df.columns:
                func_col = col
                break

        # Try to find exonic function column
        exonic_func_col = None
        for col in ['ExonicFunc.refGeneWithVer', 'ExonicFunc.refGene', 'ExonicFunc']:
            if col in self.df.columns:
                exonic_func_col = col
                break

        # Initialize categories
        categories = {
            'intergenic': 0,
            'intronic': 0,
            'exonic': 0,
            'splicing': 0,
            'UTR5': 0,
            'UTR3': 0,
            'ncRNA': 0,
        }

        # Main functional categories
        func_counts = {}
        if func_col and func_col in self.df.columns:
            func_counts = self.df[func_col].value_counts().to_dict()

        # Update categories if we have func_col data
        if func_counts:
            categories['intergenic'] = func_counts.get('intergenic', 0) + func_counts.get('intergenic_variant', 0)
            categories['intronic'] = func_counts.get('intronic', 0) + func_counts.get('intron_variant', 0)
            categories['exonic'] = func_counts.get('exonic', 0) + func_counts.get('exonic_variant', 0)
            categories['splicing'] = func_counts.get('splicing', 0) + func_counts.get('splice_site', 0)
            categories['UTR5'] = func_counts.get('UTR5', 0) + func_counts.get("5'UTR", 0)
            categories['UTR3'] = func_counts.get('UTR3', 0) + func_counts.get("3'UTR", 0)
            categories['ncRNA'] = func_counts.get('ncRNA_exonic', 0) + func_counts.get('ncRNA_intronic', 0)

        # Exonic sub-categories - get from ExonicFunc column
        exonic_counts = {}
        if exonic_func_col and exonic_func_col in self.df.columns:
            # Get all non-null exonic function values
            exonic_counts = self.df[self.df[exonic_func_col] != '.'][exonic_func_col].value_counts().to_dict()

        # Exonic sub-categories
        # ANNOVAR ExonicFunc categories:
        # - nonsynonymous SNV (missense - amino acid change)
        # - synonymous SNV (silent - no amino acid change)
        # - stopgain (nonsense - premature stop codon)
        # - stoploss (nonsense - stop codon removed)
        # - frameshift insertion/deletion
        # - nonframeshift insertion/deletion (in-frame indel)
        # - splicing (if also marked exonic)
        exonic_subcategories = {
            'nonsynonymous SNV': 0,  # Missense variants
            'synonymous SNV': 0,     # Silent variants
            'stopgain': 0,           # Nonsense - gain stop codon
            'stoploss': 0,           # Nonsense - lose stop codon
            'frameshift': 0,         # Frameshift insertions/deletions
            'nonframeshift': 0,      # In-frame insertions/deletions
            'splicing': 0,           # Splice site variants
            'unknown': 0             # Unknown/other variants
        }

        for key, count in exonic_counts.items():
            key_lower = str(key).lower().strip()
            if key_lower in ['.', '', 'nan', 'none']:
                continue
            elif 'nonsynonymous snv' in key_lower or key_lower == 'nonsynonymous_snv':
                exonic_subcategories['nonsynonymous SNV'] += count
            elif 'synonymous snv' in key_lower or key_lower == 'synonymous_snv':
                exonic_subcategories['synonymous SNV'] += count
            elif 'stopgain' in key_lower:
                exonic_subcategories['stopgain'] += count
            elif 'stoploss' in key_lower:
                exonic_subcategories['stoploss'] += count
            elif 'frameshift' in key_lower and 'nonframeshift' not in key_lower:
                exonic_subcategories['frameshift'] += count
            elif 'nonframeshift' in key_lower:
                exonic_subcategories['nonframeshift'] += count
            elif 'splicing' in key_lower:
                exonic_subcategories['splicing'] += count
            else:
                exonic_subcategories['unknown'] += count

        return {
            'categories': categories,
            'exonic_subcategories': exonic_subcategories
        }

    def get_all_metrics(self) -> Dict[str, Any]:
        """Get all metrics in one call"""
        return {
            'headline': self.get_headline_metrics(),
            'chromosome_distribution': self.get_chromosome_distribution(),
            'chromosome_distribution_normalized': self.get_chromosome_distribution(normalize=True),
            'gene_distribution': self.get_gene_distribution(),
            'gene_distribution_protein_altering': self.get_gene_distribution(protein_altering_only=True),
            'functional_impact': self.get_functional_impact()
        }
