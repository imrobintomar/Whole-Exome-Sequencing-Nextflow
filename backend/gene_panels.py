"""
Gene panel management for variant filtering
"""
import requests
from typing import List, Dict
from functools import lru_cache

class GenePanelManager:
    """Manages gene panels from PanelApp and custom sources"""

    PANELAPP_BASE_URL = "https://panelapp.genomicsengland.co.uk/api/v1"

    # ACMG Secondary Findings v3.2 (2023) - 81 genes
    ACMG_SF_GENES = [
        'BRCA1', 'BRCA2', 'MLH1', 'MSH2', 'MSH6', 'PMS2', 'APC', 'TP53',
        'MUTYH', 'PTEN', 'STK11', 'BMPR1A', 'SMAD4', 'VHL', 'MEN1', 'RET',
        'RB1', 'NF2', 'TSC1', 'TSC2', 'WT1', 'SDHD', 'SDHAF2', 'SDHC', 'SDHB',
        'FH', 'MAX', 'TMEM127', 'OTC', 'RPE65', 'BTD', 'SCN5A', 'RYR2',
        'KCNQ1', 'KCNH2', 'SCN5A', 'RYR2', 'LDLR', 'APOB', 'PCSK9', 'MYH7',
        'MYBPC3', 'TNNI3', 'TNNT2', 'TPM1', 'MYL2', 'MYL3', 'ACTC1', 'PRKAG2',
        'GLA', 'LAMP2', 'PKP2', 'DSP', 'DSC2', 'TMEM43', 'DSG2', 'RYR2',
        'SCN5A', 'KCNQ1', 'KCNH2', 'COL3A1', 'FBN1', 'TGFBR1', 'TGFBR2',
        'SMAD3', 'ACTA2', 'MYH11', 'MYLK', 'RET', 'SDHD', 'SDHAF2', 'SDHC',
        'SDHB', 'MAX', 'TMEM127', 'VHL', 'ATP7B', 'OTC', 'PCCA', 'PCCB',
        'RPE65'
    ]

    @classmethod
    @lru_cache(maxsize=100)
    def search_panels(cls, query: str) -> List[Dict]:
        """Search PanelApp for gene panels"""
        try:
            response = requests.get(
                f"{cls.PANELAPP_BASE_URL}/panels/",
                params={"search": query},
                timeout=10
            )
            response.raise_for_status()
            return response.json()['results']
        except Exception as e:
            print(f"Error fetching panels: {e}")
            return []

    @classmethod
    @lru_cache(maxsize=100)
    def get_panel_genes(cls, panel_id: int, confidence_level: str = "3") -> List[str]:
        """
        Get genes from a PanelApp panel

        Args:
            panel_id: PanelApp panel ID
            confidence_level: "3" (green/definitive), "2" (amber), "1" (red)

        Returns:
            List of gene symbols
        """
        try:
            response = requests.get(
                f"{cls.PANELAPP_BASE_URL}/panels/{panel_id}/",
                timeout=10
            )
            response.raise_for_status()
            panel_data = response.json()

            genes = []
            for gene in panel_data.get('genes', []):
                if gene['confidence_level'] == confidence_level:
                    genes.append(gene['gene_data']['gene_symbol'])

            return genes
        except Exception as e:
            print(f"Error fetching panel {panel_id}: {e}")
            return []

    @classmethod
    def get_acmg_secondary_findings_genes(cls) -> List[str]:
        """Get ACMG Secondary Findings v3.2 gene list"""
        return cls.ACMG_SF_GENES.copy()

    @classmethod
    def filter_variants_by_genes(cls, variants_df, genes: List[str], gene_column: str = 'Gene.refGeneWithVer'):
        """
        Filter variants DataFrame to only include specified genes

        Args:
            variants_df: pandas DataFrame with variants
            genes: List of gene symbols to include
            gene_column: Name of gene column in DataFrame

        Returns:
            Filtered DataFrame
        """
        import pandas as pd

        if gene_column not in variants_df.columns:
            raise ValueError(f"Column {gene_column} not found in variants")

        # Handle multiple genes per variant (separated by comma)
        def gene_in_list(gene_str):
            if pd.isna(gene_str):
                return False
            variant_genes = [g.strip() for g in str(gene_str).split(',')]
            return any(g in genes for g in variant_genes)

        mask = variants_df[gene_column].apply(gene_in_list)
        return variants_df[mask]


# Example usage functions
def filter_to_cardiac_genes(variants_df):
    """Filter to cardiac disease genes"""
    manager = GenePanelManager()
    cardiac_genes = manager.get_panel_genes(panel_id=93)  # Cardiac arrhythmias
    return manager.filter_variants_by_genes(variants_df, cardiac_genes)


def filter_to_acmg_sf(variants_df):
    """Filter to ACMG Secondary Findings genes"""
    manager = GenePanelManager()
    acmg_genes = manager.get_acmg_secondary_findings_genes()
    return manager.filter_variants_by_genes(variants_df, acmg_genes)
