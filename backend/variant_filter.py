"""
Variant filtering module for post-annotation filtering with custom parameters.

Provides deterministic, production-safe filtering of ANNOVAR-annotated variants
without re-running annotation.
"""

import pandas as pd
from pandas import DataFrame
from pydantic import BaseModel


# =========================
# Configuration
# =========================

class FilterConfig(BaseModel):
    # Functional filters
    include_exonic: bool = True
    include_splicing: bool = True
    exclude_synonymous: bool = True

    # Clinical filters
    exclude_benign: bool = True
    exclude_likely_benign: bool = True

    # Population frequency
    max_gnomad_af: float = 0.01  # UI value ignored for range, kept for API compatibility
    gnomad_af_column: str = "gnomad40_exome_AF"

    # Quality
    min_depth: int = 1
    require_pass: bool = False

    # Reserved (future)
    include_pathogenic: bool = False
    include_likely_pathogenic: bool = False
    include_vus: bool = False


# =========================
# I/O
# =========================

def load_annotated_file(file_path: str) -> DataFrame:
    """Load tab-delimited ANNOVAR output with full fidelity."""
    return pd.read_csv(
        file_path,
        sep="\t",
        keep_default_na=False,
        dtype=str,
        low_memory=False,
    )


def save_filtered_csv(df: DataFrame, output_path: str) -> None:
    """Save filtered variants as comma-delimited CSV."""
    df.to_csv(
        output_path,
        index=False,
        sep=",",
        na_rep="",
        quoting=1,  # csv.QUOTE_MINIMAL
    )


# =========================
# Normalization
# =========================

def normalize_columns(df: DataFrame) -> DataFrame:
    """
    Normalize selected annotation columns once to ensure safe comparisons.
    """
    cols = [
        "Func.refGeneWithVer",
        "ExonicFunc.refGene",
        "ExonicFunc.refGeneWithVer",
        "CLNSIG",
        "FILTER",
    ]

    for col in cols:
        if col in df.columns:
            df[col] = df[col].str.strip().str.lower()

    return df


# =========================
# Filtering Logic
# =========================

def apply_filters(df: DataFrame, config: FilterConfig) -> DataFrame:
    """
    Apply all filters sequentially.
    """

    # ---------- 1. Functional region ----------
    if config.include_exonic or config.include_splicing:
        if "Func.refGeneWithVer" in df.columns:
            patterns = []
            if config.include_exonic:
                patterns.append("exonic")
            if config.include_splicing:
                patterns.append("splicing")

            if patterns:
                df = df[
                    df["Func.refGeneWithVer"].str.contains("|".join(patterns), na=False)
                ]

    # ---------- 2. Exclude synonymous SNVs ----------
    if config.exclude_synonymous:
        exonic_col = None
        if "ExonicFunc.refGene" in df.columns:
            exonic_col = "ExonicFunc.refGene"
        elif "ExonicFunc.refGeneWithVer" in df.columns:
            exonic_col = "ExonicFunc.refGeneWithVer"

        if exonic_col:
            df = df[
                ~df[exonic_col].isin({"synonymous SNV"})
            ]

    # ---------- 3. Exclude benign / likely benign ----------
    if config.exclude_benign or config.exclude_likely_benign:
        if "CLNSIG" in df.columns:
            patterns = []
            if config.exclude_benign:
                patterns.append("benign")
            if config.exclude_likely_benign:
                patterns.append("likely benign")

            if patterns:
                mask = df["CLNSIG"].str.contains("|".join(patterns), na=False)
                df = df[~mask]

    # ---------- 4. Population frequency (0.00001 ≤ AF ≤ 0.1) ----------
    if config.gnomad_af_column in df.columns:
        MIN_AF = 0.00001
        MAX_AF = 0.1

        def af_pass(val: str) -> bool:
            try:
                af = float(val)
                return MIN_AF <= af <= MAX_AF
            except Exception:
                return val in {"", ".", "nan"}

        df = df[df[config.gnomad_af_column].apply(af_pass)]

    # ---------- 5. Read depth (DP ≥ max(1, user value)) ----------
    if "DP" in df.columns:
        effective_min_dp = max(1, config.min_depth)

        def depth_pass(val: str) -> bool:
            try:
                return float(val) >= effective_min_dp
            except Exception:
                return True

        df = df[df["DP"].apply(depth_pass)]

    # ---------- 6. PASS-only ----------
    if config.require_pass and "FILTER" in df.columns:
        df = df[df["FILTER"] == "pass"]

    return df


# =========================
# Statistics
# =========================

def get_filter_statistics(original_df: DataFrame, filtered_df: DataFrame) -> dict:
    """Generate filtering summary statistics."""
    total = len(original_df)
    kept = len(filtered_df)
    removed = total - kept
    reduction = round((removed / total) * 100, 2) if total else 0

    top_genes = []
    if "Gene.refGeneWithVer" in filtered_df.columns:
        counts = filtered_df["Gene.refGeneWithVer"].value_counts().head(10)
        top_genes = [{"gene": g, "count": int(c)} for g, c in counts.items()]

    functional_distribution = {}
    if "Func.refGeneWithVer" in filtered_df.columns:
        counts = filtered_df["Func.refGeneWithVer"].value_counts()
        functional_distribution = {k: int(v) for k, v in counts.items()}

    return {
        "total_variants": total,
        "filtered_variants": kept,
        "removed_variants": removed,
        "reduction_percentage": reduction,
        "top_genes": top_genes,
        "functional_distribution": functional_distribution,
    }
