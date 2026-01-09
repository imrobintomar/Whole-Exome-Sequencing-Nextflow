#!/usr/bin/env python3
"""
Add UniqueID column to ANNOVAR multianno.txt file
Preserves original format, just adds chr:pos:ref:alt identifier as first column
"""
import sys
import pandas as pd
from pathlib import Path

def add_unique_id(input_file, output_file):
    """
    Add UniqueID column to ANNOVAR file

    Args:
        input_file: Path to ANNOVAR .hg38_multianno.txt file
        output_file: Path to output file with UniqueID added
    """
    try:
        # Read ANNOVAR file with robust error handling
        print(f"Reading ANNOVAR file: {input_file}")

        # Try UTF-8 first, fallback to latin1 if needed
        try:
            df = pd.read_csv(
                input_file,
                sep='\t',
                encoding='utf-8',
                low_memory=False,
                on_bad_lines='skip'
            )
        except UnicodeDecodeError:
            print("UTF-8 decoding failed, trying latin1...")
            df = pd.read_csv(
                input_file,
                sep='\t',
                encoding='latin1',
                low_memory=False,
                on_bad_lines='skip'
            )

        print(f"Successfully read {len(df)} variants")
        print(f"Columns: {len(df.columns)}")

        # Validate required columns exist
        required_cols = ['Chr', 'Start', 'Ref', 'Alt']
        missing = [col for col in required_cols if col not in df.columns]
        if missing:
            raise ValueError(f"Missing required columns: {missing}")

        # Create UniqueID as chr:pos:ref:alt
        print("Creating UniqueID column...")
        df.insert(0, 'UniqueID',
                  df['Chr'].astype(str) + ':' +
                  df['Start'].astype(str) + ':' +
                  df['Ref'].astype(str) + ':' +
                  df['Alt'].astype(str))

        print(f"Created UniqueID for {len(df)} variants")

        # Write to output file with proper encoding (keeping tab-delimited format)
        df.to_csv(
            output_file,
            sep='\t',
            index=False,
            encoding='utf-8'
        )

        print(f"Successfully wrote {len(df)} variants to {output_file}")
        print(f"Output file size: {Path(output_file).stat().st_size / 1024:.2f} KB")

        return 0

    except Exception as e:
        print(f"ERROR: Failed to add UniqueID: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: add_unique_id.py <input.txt> <output.txt>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    sys.exit(add_unique_id(input_file, output_file))
