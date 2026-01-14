"""
Variant filtering module for post-annotation filtering with custom parameters.

This module provides deterministic, production-safe filtering of
ANNOVAR-annotated variant files without re-running annotation.
"""

import pandas as pd
from pandas import DataFrame
from pydantic import BaseModel


# ============================================================
# Configuration Model
# ============================================================

class FilterConfig(BaseModel):
    # Functional filters
    include_exonic: bool = True
    include_splicing: bool = True
    exclude_synonymous: bool = True

    # Clinical significance filters
    exclude_benign: bool = True
    exclude_likely_benign: bool = True

    # Population frequency filters
    # Only upper bound is enforced (keep rare variants)
    max_gnomad_af: float = 0.1
    gnomad_af_column: str = "gnomad40_exome_AF"

    # Quality filters
    min_depth: int = 1
    require_pass: bool = False

    # Reserved for future ACMG logic
    include_pathogenic: bool = False
    include_likely_pathogenic: bool = False
    include_vus: bool = False


# ============================================================
# I/O Utilities
# ============================================================

def load_annotated_file(file_path: str) -> DataFrame:
    """
    Load ANNOVAR annotated file (tab-delimited) with full fidelity.
    """
    return pd.read_csv(
        file_path,
        sep="\t",
        keep_default_na=False,
        dtype=str,
        low_memory=False
    )


def save_filtered_csv(df: DataFrame, output_path: str) -> None:
    """
    Save filtered DataFrame as comma-delimited CSV.
    """
    df.to_csv(
        output_path,
        index=False,
        sep=",",
        na_rep="",
        quoting=1  # csv.QUOTE_MINIMAL
    )


# ============================================================
# Normalization
# ============================================================

def normalize_columns(df: DataFrame) -> DataFrame:
    """
    Normalize selected columns once to ensure safe string comparisons.
    """
    cols = [
        "Func.refGeneWithVer",
        "ExonicFunc.refGeneWithVer",
        "ExonicFunc.refGene",
        "CLNSIG",
        "FILTER",
    ]

    for col in cols:
        if col in df.columns:
            df[col] = df[col].str.strip().str.lower()

    return df


# ============================================================
# Core Filtering Logic
# ============================================================

def apply_filters(df: DataFrame, config: FilterConfig) -> DataFrame:
    """
    Apply all enabled filters sequentially.
    """

    # --------------------------------------------------------
    # 1. Functional region (exonic / splicing)
    # --------------------------------------------------------
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

    # --------------------------------------------------------
    # 2. Exclude synonymous SNVs ONLY
    # --------------------------------------------------------
    if config.exclude_synonymous:
        exonic_col = None
        if "ExonicFunc.refGeneWithVer" in df.columns:
            exonic_col = "ExonicFunc.refGeneWithVer"
        elif "ExonicFunc.refGene" in df.columns:
            exonic_col = "ExonicFunc.refGene"

        if exonic_col:
            df = df[
                ~df[exonic_col].isin({"synonymous snv"})
            ]

    # --------------------------------------------------------
    # 3. Exclude benign / likely benign (ClinVar)
    # --------------------------------------------------------
    if config.exclude_benign or config.exclude_likely_benign:
        if "CLNSIG" in df.columns:
            patterns = []
            if config.exclude_benign:
                patterns.append("benign")
            if config.exclude_likely_benign:
                patterns.append("likely benign")

            if patterns:
                df = df[
                    ~df["CLNSIG"].str.contains("|".join(patterns), na=False)
                ]

    # --------------------------------------------------------
    # 4. Population frequency (KEEP ALL RARE VARIANTS)
    #     Rule: AF ≤ max_gnomad_af
    #     NO LOWER BOUND
    # --------------------------------------------------------
    if config.gnomad_af_column in df.columns:
        MAX_AF = config.max_gnomad_af

        def af_pass(val: str) -> bool:
            try:
                return float(val) <= MAX_AF
            except Exception:
                # Missing AF (., empty) → keep
                return True

        df = df[df[config.gnomad_af_column].apply(af_pass)]

    # --------------------------------------------------------
    # 5. Read depth (DP ≥ max(1, user value))
    # --------------------------------------------------------
    if "DP" in df.columns:
        effective_min_dp = max(1, config.min_depth)

        def depth_pass(val: str) -> bool:
            try:
                return float(val) >= effective_min_dp
            except Exception:
                return True  # keep if missing

        df = df[df["DP"].apply(depth_pass)]

    # --------------------------------------------------------
    # 6. PASS-only filter
    # --------------------------------------------------------
    if config.require_pass and "FILTER" in df.columns:
        df = df[df["FILTER"] == "pass"]

    return df


# ============================================================
# Statistics
# ============================================================

def get_filter_statistics(original_df: DataFrame, filtered_df: DataFrame) -> dict:
    """
    Compute summary statistics comparing original and filtered datasets.
    """
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
