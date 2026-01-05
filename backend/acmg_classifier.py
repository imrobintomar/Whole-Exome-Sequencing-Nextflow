"""
ACMG/AMP 2015 Guidelines Variant Classification Engine

Implements automated classification of variants according to:
Richards et al. (2015) Standards and guidelines for the interpretation of
sequence variants: a joint consensus recommendation of the American College
of Medical Genetics and Genomics and the Association for Molecular Pathology.

Classification: Pathogenic (P), Likely Pathogenic (LP), VUS, Likely Benign (LB), Benign (B)
"""

from enum import Enum
from typing import Dict, List, Optional, Set
from dataclasses import dataclass
import re


class ACMGStrength(Enum):
    """Evidence strength levels"""
    VERY_STRONG = "very_strong"
    STRONG = "strong"
    MODERATE = "moderate"
    SUPPORTING = "supporting"
    STAND_ALONE = "stand_alone"


class ACMGClassification(Enum):
    """Final ACMG classifications"""
    PATHOGENIC = "Pathogenic"
    LIKELY_PATHOGENIC = "Likely Pathogenic"
    VUS = "Uncertain Significance"
    LIKELY_BENIGN = "Likely Benign"
    BENIGN = "Benign"


@dataclass
class ACMGEvidence:
    """Evidence for a specific ACMG criterion"""
    code: str
    met: bool
    strength: ACMGStrength
    description: str
    auto_applied: bool = True  # Was it auto-applied or manually curated?


class ACMGClassifier:
    """
    ACMG/AMP variant classification engine

    Pathogenic criteria:
    - PVS1: Null variant in LOF-intolerant gene
    - PS1-PS4: Strong evidence
    - PM1-PM6: Moderate evidence
    - PP1-PP5: Supporting evidence

    Benign criteria:
    - BA1: Stand-alone benign
    - BS1-BS4: Strong benign
    - BP1-BP7: Supporting benign
    """

    # Classification combinations (from ACMG 2015 guidelines)
    PATHOGENIC_COMBINATIONS = [
        # Very strong + other
        {"PVS1": 1, "PS": 1},
        {"PVS1": 1, "PM": 2},
        {"PVS1": 1, "PM": 1, "PP": 2},
        {"PVS1": 1, "PP": 4},
        # Strong only
        {"PS": 2},
        {"PS": 1, "PM": 3},
        {"PS": 1, "PM": 2, "PP": 2},
        {"PS": 1, "PM": 1, "PP": 4},
    ]

    LIKELY_PATHOGENIC_COMBINATIONS = [
        {"PVS1": 1, "PM": 1},
        {"PS": 1, "PM": 1},
        {"PS": 1, "PM": 2},
        {"PS": 1, "PP": 2},
        {"PM": 3},
        {"PM": 2, "PP": 2},
        {"PM": 1, "PP": 4},
    ]

    BENIGN_COMBINATIONS = [
        {"BA1": 1},  # Stand-alone
    ]

    LIKELY_BENIGN_COMBINATIONS = [
        {"BS": 1, "BP": 1},
        {"BP": 2},
    ]

    def __init__(self):
        self.evidence: Dict[str, ACMGEvidence] = {}

    def classify_variant(self, variant: Dict) -> Dict:
        """
        Main classification function

        Args:
            variant: Dictionary with variant information
                Required fields:
                - consequence: Effect (missense, nonsense, frameshift, etc.)
                - gene: Gene symbol
                - af_gnomad: gnomAD allele frequency
                - cadd_phred: CADD score
                - revel_score: REVEL score
                - sift_pred: SIFT prediction
                - polyphen_pred: PolyPhen prediction
                - clinvar_sig: ClinVar significance (if available)

                Optional fields:
                - pli: Gene pLI score
                - loeuf: Gene LOEUF score
                - spliceai_max: Max SpliceAI score
                - in_protein_domain: Boolean

        Returns:
            Dictionary with classification and evidence
        """
        self.evidence = {}

        # Apply pathogenic criteria
        self._apply_pvs1(variant)
        self._apply_ps1(variant)
        self._apply_pm1(variant)
        self._apply_pm2(variant)
        self._apply_pm5(variant)
        self._apply_pp2(variant)
        self._apply_pp3(variant)
        self._apply_pp5(variant)

        # Apply benign criteria
        self._apply_ba1(variant)
        self._apply_bs1(variant)
        self._apply_bp1(variant)
        self._apply_bp3(variant)
        self._apply_bp4(variant)
        self._apply_bp7(variant)

        # Calculate final classification
        classification = self._determine_classification()

        # Count evidence by strength
        evidence_summary = self._summarize_evidence()

        return {
            "classification": classification.value,
            "evidence": {code: {
                "met": ev.met,
                "strength": ev.strength.value,
                "description": ev.description,
                "auto_applied": ev.auto_applied
            } for code, ev in self.evidence.items() if ev.met},
            "evidence_summary": evidence_summary,
            "met_criteria": [code for code, ev in self.evidence.items() if ev.met]
        }

    # ========== PATHOGENIC CRITERIA ==========

    def _apply_pvs1(self, variant: Dict):
        """
        PVS1: Null variant (nonsense, frameshift, splice) in gene where LOF is disease mechanism
        """
        consequence = variant.get("consequence", "").lower()
        pli = variant.get("pli", 0)
        loeuf = variant.get("loeuf", 1.0)

        # Check if null variant
        is_null = any(term in consequence for term in [
            "stop_gained", "nonsense", "frameshift",
            "splice_donor", "splice_acceptor"
        ])

        # Check if gene is LOF-intolerant (pLI > 0.9 or LOEUF < 0.35)
        is_lof_intolerant = pli > 0.9 or loeuf < 0.35

        met = is_null and is_lof_intolerant

        description = ""
        if met:
            description = f"Null variant ({consequence}) in LOF-intolerant gene (pLI={pli:.2f}, LOEUF={loeuf:.2f})"
        elif is_null and not is_lof_intolerant:
            description = f"Null variant but gene is LOF-tolerant (pLI={pli:.2f})"

        self.evidence["PVS1"] = ACMGEvidence(
            code="PVS1",
            met=met,
            strength=ACMGStrength.VERY_STRONG,
            description=description
        )

    def _apply_ps1(self, variant: Dict):
        """
        PS1: Same amino acid change as established pathogenic variant
        """
        # This requires database of known pathogenic variants
        # Simplified: Check if in ClinVar as pathogenic
        clinvar = variant.get("clinvar_sig", "").lower()

        met = "pathogenic" in clinvar and "conflicting" not in clinvar

        self.evidence["PS1"] = ACMGEvidence(
            code="PS1",
            met=met,
            strength=ACMGStrength.STRONG,
            description=f"ClinVar: {variant.get('clinvar_sig', 'N/A')}" if met else ""
        )

    def _apply_pm1(self, variant: Dict):
        """
        PM1: Located in mutational hot spot or critical functional domain
        """
        in_domain = variant.get("in_protein_domain", False)

        self.evidence["PM1"] = ACMGEvidence(
            code="PM1",
            met=in_domain,
            strength=ACMGStrength.MODERATE,
            description="Variant in critical protein domain" if in_domain else ""
        )

    def _apply_pm2(self, variant: Dict):
        """
        PM2: Absent or extremely low frequency in population databases (gnomAD)

        Thresholds:
        - Recessive: AF < 0.01 (1%)
        - Dominant: AF < 0.0001 (0.01%)
        """
        af = variant.get("af_gnomad", 0)

        # Using dominant threshold (more stringent)
        met = af < 0.0001

        self.evidence["PM2"] = ACMGEvidence(
            code="PM2",
            met=met,
            strength=ACMGStrength.MODERATE,
            description=f"gnomAD AF = {af:.6f} (< 0.0001 threshold)" if met else f"gnomAD AF = {af:.6f}"
        )

    def _apply_pm5(self, variant: Dict):
        """
        PM5: Novel missense at amino acid where different pathogenic missense exists
        """
        # Requires database lookup - simplified here
        # Could check ClinVar for other missense at same position

        self.evidence["PM5"] = ACMGEvidence(
            code="PM5",
            met=False,
            strength=ACMGStrength.MODERATE,
            description="Requires database of known variants"
        )

    def _apply_pp2(self, variant: Dict):
        """
        PP2: Missense in gene with low rate of benign missense variation
        """
        consequence = variant.get("consequence", "").lower()
        loeuf = variant.get("loeuf", 1.0)

        is_missense = "missense" in consequence
        low_missense_rate = loeuf < 0.5  # Gene intolerant to variation

        met = is_missense and low_missense_rate

        self.evidence["PP2"] = ACMGEvidence(
            code="PP2",
            met=met,
            strength=ACMGStrength.SUPPORTING,
            description=f"Missense in constrained gene (LOEUF={loeuf:.2f})" if met else ""
        )

    def _apply_pp3(self, variant: Dict):
        """
        PP3: Multiple computational evidence support deleteriousness

        Criteria: ≥4 of 5 predictors damaging:
        - CADD > 20
        - REVEL > 0.5
        - SIFT = Deleterious
        - PolyPhen = Damaging
        - SpliceAI > 0.5 (for splice variants)
        """
        cadd = variant.get("cadd_phred", 0)
        revel = variant.get("revel_score", 0)
        sift = variant.get("sift_pred", "").lower()
        polyphen = variant.get("polyphen_pred", "").lower()
        spliceai = variant.get("spliceai_max", 0)

        damaging_count = 0
        predictors = []

        if cadd > 20:
            damaging_count += 1
            predictors.append(f"CADD={cadd:.1f}")
        if revel > 0.5:
            damaging_count += 1
            predictors.append(f"REVEL={revel:.3f}")
        if "deleterious" in sift or "damaging" in sift:
            damaging_count += 1
            predictors.append(f"SIFT={sift}")
        if "damaging" in polyphen or "probably_damaging" in polyphen:
            damaging_count += 1
            predictors.append(f"PolyPhen={polyphen}")
        if spliceai > 0.5:
            damaging_count += 1
            predictors.append(f"SpliceAI={spliceai:.3f}")

        met = damaging_count >= 4

        self.evidence["PP3"] = ACMGEvidence(
            code="PP3",
            met=met,
            strength=ACMGStrength.SUPPORTING,
            description=f"{damaging_count}/5 predictors damaging: {', '.join(predictors)}" if predictors else ""
        )

    def _apply_pp5(self, variant: Dict):
        """
        PP5: Reputable source reports variant as pathogenic
        """
        clinvar = variant.get("clinvar_sig", "").lower()

        # Check for pathogenic assertions (but not conflicting)
        met = "pathogenic" in clinvar and "conflicting" not in clinvar

        # Downgrade from PS1 if already applied
        if met and self.evidence.get("PS1", ACMGEvidence("", False, ACMGStrength.SUPPORTING, "")).met:
            met = False  # Don't double-count

        self.evidence["PP5"] = ACMGEvidence(
            code="PP5",
            met=met,
            strength=ACMGStrength.SUPPORTING,
            description=f"ClinVar: {variant.get('clinvar_sig', '')}" if met else ""
        )

    # ========== BENIGN CRITERIA ==========

    def _apply_ba1(self, variant: Dict):
        """
        BA1: Stand-alone benign - Allele frequency > 5% in population
        """
        af = variant.get("af_gnomad", 0)

        met = af > 0.05

        self.evidence["BA1"] = ACMGEvidence(
            code="BA1",
            met=met,
            strength=ACMGStrength.STAND_ALONE,
            description=f"gnomAD AF = {af:.4f} (>5%)" if met else ""
        )

    def _apply_bs1(self, variant: Dict):
        """
        BS1: Allele frequency greater than expected for disorder

        Simplified: AF > 1% for recessive, > 0.1% for dominant
        """
        af = variant.get("af_gnomad", 0)

        # Using recessive threshold
        met = af > 0.01 and af <= 0.05  # BA1 takes precedence if >5%

        self.evidence["BS1"] = ACMGEvidence(
            code="BS1",
            met=met,
            strength=ACMGStrength.STRONG,
            description=f"gnomAD AF = {af:.4f} (>1%)" if met else ""
        )

    def _apply_bp1(self, variant: Dict):
        """
        BP1: Missense in gene where only truncating variants cause disease
        """
        consequence = variant.get("consequence", "").lower()
        pli = variant.get("pli", 0)

        is_missense = "missense" in consequence
        lof_mechanism = pli > 0.9  # Gene has LOF mechanism

        met = is_missense and lof_mechanism

        self.evidence["BP1"] = ACMGEvidence(
            code="BP1",
            met=met,
            strength=ACMGStrength.SUPPORTING,
            description=f"Missense in LOF gene (pLI={pli:.2f})" if met else ""
        )

    def _apply_bp3(self, variant: Dict):
        """
        BP3: In-frame indels in non-repeat region without known function
        """
        consequence = variant.get("consequence", "").lower()

        is_inframe = "inframe" in consequence

        # Simplified - would need repeat region annotation
        met = is_inframe

        self.evidence["BP3"] = ACMGEvidence(
            code="BP3",
            met=met,
            strength=ACMGStrength.SUPPORTING,
            description="In-frame indel" if met else ""
        )

    def _apply_bp4(self, variant: Dict):
        """
        BP4: Multiple computational evidence support benign

        Opposite of PP3: ≥4 of 5 predictors benign
        """
        cadd = variant.get("cadd_phred", 0)
        revel = variant.get("revel_score", 0)
        sift = variant.get("sift_pred", "").lower()
        polyphen = variant.get("polyphen_pred", "").lower()

        benign_count = 0
        predictors = []

        if cadd < 15:
            benign_count += 1
            predictors.append(f"CADD={cadd:.1f}")
        if revel < 0.3:
            benign_count += 1
            predictors.append(f"REVEL={revel:.3f}")
        if "tolerated" in sift:
            benign_count += 1
            predictors.append(f"SIFT={sift}")
        if "benign" in polyphen:
            benign_count += 1
            predictors.append(f"PolyPhen={polyphen}")

        met = benign_count >= 3

        self.evidence["BP4"] = ACMGEvidence(
            code="BP4",
            met=met,
            strength=ACMGStrength.SUPPORTING,
            description=f"{benign_count}/4 predictors benign: {', '.join(predictors)}" if predictors else ""
        )

    def _apply_bp7(self, variant: Dict):
        """
        BP7: Synonymous variant with no predicted splice impact
        """
        consequence = variant.get("consequence", "").lower()
        spliceai = variant.get("spliceai_max", 0)

        is_synonymous = "synonymous" in consequence
        no_splice_impact = spliceai < 0.2

        met = is_synonymous and no_splice_impact

        self.evidence["BP7"] = ACMGEvidence(
            code="BP7",
            met=met,
            strength=ACMGStrength.SUPPORTING,
            description=f"Synonymous, SpliceAI={spliceai:.3f}" if met else ""
        )

    # ========== CLASSIFICATION LOGIC ==========

    def _summarize_evidence(self) -> Dict[str, int]:
        """Count evidence by category"""
        summary = {
            "PVS": 0,
            "PS": 0,
            "PM": 0,
            "PP": 0,
            "BA": 0,
            "BS": 0,
            "BP": 0
        }

        for code, evidence in self.evidence.items():
            if evidence.met:
                if code.startswith("PVS"):
                    summary["PVS"] += 1
                elif code.startswith("PS"):
                    summary["PS"] += 1
                elif code.startswith("PM"):
                    summary["PM"] += 1
                elif code.startswith("PP"):
                    summary["PP"] += 1
                elif code.startswith("BA"):
                    summary["BA"] += 1
                elif code.startswith("BS"):
                    summary["BS"] += 1
                elif code.startswith("BP"):
                    summary["BP"] += 1

        return summary

    def _determine_classification(self) -> ACMGClassification:
        """
        Determine final classification based on met criteria
        """
        summary = self._summarize_evidence()

        # Check stand-alone benign first
        if summary["BA"] >= 1:
            return ACMGClassification.BENIGN

        # Check pathogenic combinations
        if self._meets_pathogenic_criteria(summary):
            return ACMGClassification.PATHOGENIC

        # Check likely pathogenic
        if self._meets_likely_pathogenic_criteria(summary):
            return ACMGClassification.LIKELY_PATHOGENIC

        # Check benign
        if self._meets_benign_criteria(summary):
            return ACMGClassification.BENIGN

        # Check likely benign
        if self._meets_likely_benign_criteria(summary):
            return ACMGClassification.LIKELY_BENIGN

        # Default to VUS
        return ACMGClassification.VUS

    def _meets_pathogenic_criteria(self, summary: Dict[str, int]) -> bool:
        """Check if meets pathogenic classification"""
        # PVS1 + 1 PS
        if summary["PVS"] >= 1 and summary["PS"] >= 1:
            return True
        # PVS1 + 2 PM
        if summary["PVS"] >= 1 and summary["PM"] >= 2:
            return True
        # PVS1 + 1 PM + 2 PP
        if summary["PVS"] >= 1 and summary["PM"] >= 1 and summary["PP"] >= 2:
            return True
        # PVS1 + 4 PP
        if summary["PVS"] >= 1 and summary["PP"] >= 4:
            return True
        # 2 PS
        if summary["PS"] >= 2:
            return True
        # 1 PS + 3 PM
        if summary["PS"] >= 1 and summary["PM"] >= 3:
            return True
        # 1 PS + 2 PM + 2 PP
        if summary["PS"] >= 1 and summary["PM"] >= 2 and summary["PP"] >= 2:
            return True
        # 1 PS + 1 PM + 4 PP
        if summary["PS"] >= 1 and summary["PM"] >= 1 and summary["PP"] >= 4:
            return True

        return False

    def _meets_likely_pathogenic_criteria(self, summary: Dict[str, int]) -> bool:
        """Check if meets likely pathogenic classification"""
        # PVS1 + 1 PM
        if summary["PVS"] >= 1 and summary["PM"] >= 1:
            return True
        # 1 PS + 1-2 PM
        if summary["PS"] >= 1 and summary["PM"] >= 1:
            return True
        # 1 PS + 2 PP
        if summary["PS"] >= 1 and summary["PP"] >= 2:
            return True
        # 3 PM
        if summary["PM"] >= 3:
            return True
        # 2 PM + 2 PP
        if summary["PM"] >= 2 and summary["PP"] >= 2:
            return True
        # 1 PM + 4 PP
        if summary["PM"] >= 1 and summary["PP"] >= 4:
            return True

        return False

    def _meets_benign_criteria(self, summary: Dict[str, int]) -> bool:
        """Check if meets benign classification"""
        # 2 BS
        if summary["BS"] >= 2:
            return True

        return False

    def _meets_likely_benign_criteria(self, summary: Dict[str, int]) -> bool:
        """Check if meets likely benign classification"""
        # 1 BS + 1 BP
        if summary["BS"] >= 1 and summary["BP"] >= 1:
            return True
        # 2 BP
        if summary["BP"] >= 2:
            return True

        return False


# Example usage
if __name__ == "__main__":
    # Test case: Pathogenic variant
    test_variant = {
        "consequence": "stop_gained",
        "gene": "BRCA1",
        "af_gnomad": 0.000001,
        "cadd_phred": 35,
        "revel_score": 0.9,
        "sift_pred": "deleterious",
        "polyphen_pred": "probably_damaging",
        "clinvar_sig": "Pathogenic",
        "pli": 0.98,
        "loeuf": 0.25,
    }

    classifier = ACMGClassifier()
    result = classifier.classify_variant(test_variant)

    print(f"Classification: {result['classification']}")
    print(f"Evidence summary: {result['evidence_summary']}")
    print(f"Met criteria: {', '.join(result['met_criteria'])}")
    for code, evidence in result['evidence'].items():
        print(f"  {code}: {evidence['description']}")
