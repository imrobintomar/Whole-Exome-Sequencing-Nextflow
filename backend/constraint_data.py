"""
Gene constraint data loader for ACMG classification

Loads gnomAD gene constraint metrics (pLI, LOEUF) for determining
gene loss-of-function intolerance (used in PVS1, PP2, BP1 criteria)
"""

import pandas as pd
from pathlib import Path
from typing import Dict, Optional
from functools import lru_cache


class GeneConstraintDB:
    """Loads and queries gene constraint metrics"""

    def __init__(self, constraint_file: Optional[str] = None):
        """
        Initialize constraint database

        Args:
            constraint_file: Path to gnomAD constraint TSV file
                           If None, uses default location: data/gnomad_constraint.tsv
        """
        if constraint_file is None:
            constraint_file = Path(__file__).parent / "data" / "gnomad_constraint.tsv"

        self.constraint_file = Path(constraint_file)
        self._data: Optional[pd.DataFrame] = None
        self._gene_index: Optional[Dict] = None

    def load(self) -> bool:
        """
        Load constraint data from file

        Returns:
            True if loaded successfully, False otherwise
        """
        if not self.constraint_file.exists():
            print(f" Constraint file not found: {self.constraint_file}")
            print("Run: bash backend/download_gnomad_constraint.sh")
            return False

        try:
            # Load TSV file
            self._data = pd.read_csv(self.constraint_file, sep='\t')

            # Create gene symbol index for fast lookup
            if 'gene' in self._data.columns:
                self._gene_index = {}
                for idx, row in self._data.iterrows():
                    gene = row['gene']
                    self._gene_index[gene] = {
                        'pli': row.get('pLI', row.get('oe_lof_upper', 0)),
                        'loeuf': row.get('oe_lof_upper', row.get('LOEUF', 1.0)),
                        'oe_mis': row.get('oe_mis_upper', 1.0),
                        'syn_z': row.get('syn_z', 0)
                    }

                print(f" Loaded constraint data for {len(self._gene_index)} genes")
                return True
            else:
                print(f" 'gene' column not found in {self.constraint_file}")
                return False

        except Exception as e:
            print(f" Error loading constraint data: {e}")
            return False

    @lru_cache(maxsize=10000)
    def get_gene_constraint(self, gene: str) -> Dict[str, float]:
        """
        Get constraint metrics for a gene

        Args:
            gene: Gene symbol (e.g., "BRCA1")

        Returns:
            Dictionary with constraint metrics:
            - pli: Probability of LOF intolerance (0-1, >0.9 = intolerant)
            - loeuf: LOEUF score (<0.35 = intolerant)
            - oe_mis: Missense constraint
            - syn_z: Synonymous Z-score
        """
        if self._gene_index is None:
            if not self.load():
                return {"pli": 0, "loeuf": 1.0, "oe_mis": 1.0, "syn_z": 0}

        return self._gene_index.get(gene, {
            "pli": 0,
            "loeuf": 1.0,
            "oe_mis": 1.0,
            "syn_z": 0
        })

    def is_lof_intolerant(self, gene: str, pli_threshold: float = 0.9,
                          loeuf_threshold: float = 0.35) -> bool:
        """
        Check if gene is loss-of-function intolerant

        Args:
            gene: Gene symbol
            pli_threshold: pLI threshold (default 0.9)
            loeuf_threshold: LOEUF threshold (default 0.35)

        Returns:
            True if gene is LOF-intolerant
        """
        metrics = self.get_gene_constraint(gene)
        return metrics["pli"] > pli_threshold or metrics["loeuf"] < loeuf_threshold

    def is_missense_constrained(self, gene: str, oe_mis_threshold: float = 0.8) -> bool:
        """
        Check if gene is constrained for missense variation

        Args:
            gene: Gene symbol
            oe_mis_threshold: Observed/Expected missense threshold

        Returns:
            True if gene has low missense tolerance
        """
        metrics = self.get_gene_constraint(gene)
        return metrics["oe_mis"] < oe_mis_threshold


# Global instance
_constraint_db = None


def get_constraint_db() -> GeneConstraintDB:
    """Get or create global constraint database instance"""
    global _constraint_db
    if _constraint_db is None:
        _constraint_db = GeneConstraintDB()
        _constraint_db.load()
    return _constraint_db


# Example usage
if __name__ == "__main__":
    db = GeneConstraintDB()
    if db.load():
        # Test some well-known genes
        test_genes = ["BRCA1", "TP53", "TTN", "OR4F5"]

        for gene in test_genes:
            metrics = db.get_gene_constraint(gene)
            lof_intolerant = db.is_lof_intolerant(gene)

            print(f"\n{gene}:")
            print(f"  pLI: {metrics['pli']:.3f}")
            print(f"  LOEUF: {metrics['loeuf']:.3f}")
            print(f"  LOF intolerant: {lof_intolerant}")
