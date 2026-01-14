"""
Variant filtering module for post-annotation filtering with custom parameters.

This module provides on-demand filtering of ANNOVAR-annotated variant files
without re-running the expensive annotation pipeline.
"""

import pandas as pd
from pandas import DataFrame
from pydantic import BaseModel
from typing import Optional


class FilterConfig(BaseModel):
    """
    Configuration for variant filtering parameters.

    Users can customize these parameters to filter variants based on:
    - Functional impact (exonic/splicing)
    - Clinical significance (benign/pathogenic)
    - Population frequency (gnomAD AF)
    - Quality metrics (read depth, PASS filter)
    """

    # Functional filters
    include_exonic: bool = True
    include_splicing: bool = True
    exclude_synonymous: bool = True

    # Clinical significance filters
    exclude_benign: bool = True
    exclude_likely_benign: bool = True

    # Population frequency filters
    max_gnomad_af: float = 0.05
    gnomad_af_column: str = "gnomad40_exome_AF"

    # Quality filters
    min_depth: int = 5
    require_pass: bool = False

    # ACMG classification filters (optional)
    include_pathogenic: bool = False
    include_likely_pathogenic: bool = False
    include_vus: bool = False


def load_annotated_file(file_path: str) -> DataFrame:
    """
    Load ANNOVAR annotated file (tab-delimited Final_.txt).

    CRITICAL: Uses keep_default_na=False to preserve dots (.) and empty strings
    exactly as they appear in the original file. All values are kept as strings
    to prevent any data loss or transformation.

    Args:
        file_path: Path to the tab-delimited annotated file

    Returns:
        DataFrame with all original values preserved
    """
    return pd.read_csv(
        file_path,
        sep='\t',  # Input is tab-delimited
        keep_default_na=False,  # Don't convert '.' or '' to NaN
        dtype=str,  # Keep all values as strings to prevent data loss
        low_memory=False  # Handle large files
    )


def save_filtered_csv(df: DataFrame, output_path: str) -> None:
    """
    Save filtered DataFrame as CSV (comma-delimited) matching sample_filter-file.csv format.

    CRITICAL:
    - Changes delimiter from tab to comma
    - Preserves ALL data values exactly (dots, spaces, empty strings)
    - Doesn't convert anything to NaN
    - Quotes fields containing commas to prevent parsing issues

    Args:
        df: Filtered DataFrame to save
        output_path: Path where CSV file will be written
    """
    df.to_csv(
        output_path,
        index=False,
        sep=',',  # Output is comma-delimited
        na_rep='',  # Keep empty values as empty strings
        quoting=1  # Quote fields containing commas (csv.QUOTE_MINIMAL)
    )


def apply_filters(df: DataFrame, config: FilterConfig) -> DataFrame:
    """
    Apply all enabled filters sequentially to the variant DataFrame.

    This function filters ROWS only - it never modifies data values.
    All original values (dots, spaces, empty strings) are preserved.

    Args:
        df: DataFrame containing annotated variants
        config: FilterConfig with user-specified filter parameters

    Returns:
        Filtered DataFrame with subset of rows matching filter criteria
    """

    # Filter 1: Functional impact (exonic/splicing)
    if config.include_exonic or config.include_splicing:
        if "Func.refGeneWithVer" in df.columns:
            patterns = []
            if config.include_exonic:
                patterns.append("exonic")
            if config.include_splicing:
                patterns.append("splicing")

            if patterns:
                pattern = "|".join(patterns)
                df = df[df["Func.refGeneWithVer"].str.contains(pattern, case=False, na=False)]

    # Filter 2: Exclude synonymous SNVs
    if config.exclude_synonymous:
          exonic_col = None
          if "ExonicFunc.refGene" in df.columns:
            exonic_col = "ExonicFunc.refGene"
          elif "ExonicFunc.refGeneWithVer" in df.columns:
            exonic_col = "ExonicFunc.refGeneWithVer"
          if exonic_col:
           df = df[
             ~df[exonic_col]
             .str.strip()
             .str.lower()
             .isin({"synonymous snv"})
        ]
    # Filter 3: Exclude benign variants
    if config.exclude_benign or config.exclude_likely_benign:
        if "CLNSIG" in df.columns:
            patterns = []
            if config.exclude_benign:
                patterns.append("benign")
            if config.exclude_likely_benign:
                patterns.append("likely benign")

            if patterns:
                pattern = "|".join(patterns)
                # Create a mask for rows to exclude
                mask = df["CLNSIG"].str.lower().str.contains(pattern, na=False)
                df = df[~mask]

    # Filter 4: gnomAD allele frequency
    if config.gnomad_af_column in df.columns:
        def filter_af(val):
            """
            Filter based on allele frequency.
            Keep variants with AF < threshold.
            Treat missing values (., nan, empty) as acceptable (keep them).
            """
            try:
                af_value = float(str(val).strip())
                return af_value < config.max_gnomad_af
            except (ValueError, TypeError):
                # If conversion fails, check if it's a missing value indicator
                val_str = str(val).strip().lower()
                return val_str in ["", ".", "nan", "0", "0.0"]

        df = df[df[config.gnomad_af_column].apply(filter_af)]

    # Filter 5: Read depth (DP)
    if 'DP' in df.columns:
        # Convert DP to numeric, keeping original strings for non-numeric values
        # Then filter based on numeric comparison
        def filter_depth(val):
            try:
                dp_value = float(str(val).strip())
                return dp_value > config.min_depth
            except (ValueError, TypeError):
                # Keep rows with missing DP values
                return True

        df = df[df["DP"].apply(filter_depth)]

    # Filter 6: PASS-only variants
    if config.require_pass and 'FILTER' in df.columns:
        df = df[df["FILTER"] == "PASS"]

    return df


def get_filter_statistics(original_df: DataFrame, filtered_df: DataFrame) -> dict:
    """
    Calculate statistics comparing original and filtered variant sets.

    Args:
        original_df: Original unfiltered DataFrame
        filtered_df: Filtered DataFrame

    Returns:
        Dictionary with statistics about filtering results
    """
    total_variants = len(original_df)
    filtered_variants = len(filtered_df)
    reduction_percentage = ((total_variants - filtered_variants) / total_variants * 100) if total_variants > 0 else 0

    # Get top genes in filtered set
    top_genes = []
    if "Gene.refGeneWithVer" in filtered_df.columns:
        gene_counts = filtered_df["Gene.refGeneWithVer"].value_counts().head(10)
        top_genes = [{"gene": gene, "count": int(count)} for gene, count in gene_counts.items()]

    # Get functional distribution
    functional_distribution = {}
    if "Func.refGeneWithVer" in filtered_df.columns:
        func_counts = filtered_df["Func.refGeneWithVer"].value_counts()
        functional_distribution = {func: int(count) for func, count in func_counts.items()}

    return {
        "total_variants": total_variants,
        "filtered_variants": filtered_variants,
        "removed_variants": total_variants - filtered_variants,
        "reduction_percentage": round(reduction_percentage, 2),
        "top_genes": top_genes,
        "functional_distribution": functional_distribution
    }
