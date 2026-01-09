process addUniqueID {

    tag "$sample_id"

    input:
        tuple val(sample_id), path(annovar_txt)

    output:
        path("${sample_id}_Final_.txt")

    script:
        """
        # Validate input file exists
        if [ ! -f "${annovar_txt}" ]; then
            echo "ERROR: Input ANNOVAR file not found: ${annovar_txt}"
            exit 1
        fi

        # Add UniqueID column to ANNOVAR output
        python3 ${projectDir}/processes/add_unique_id.py \\
            "${annovar_txt}" \\
            "${sample_id}_Final_.txt"

        # Check exit code
        if [ \$? -ne 0 ]; then
            echo "ERROR: Python script failed!"
            exit 1
        fi

        # Validate output
        if [ ! -s "${sample_id}_Final_.txt" ]; then
            echo "ERROR: Output file is empty!"
            exit 1
        fi

        # Verify we have at least a header line
        LINE_COUNT=\$(wc -l < "${sample_id}_Final_.txt")
        if [ \$LINE_COUNT -lt 2 ]; then
            echo "ERROR: Output file has no data rows (only \$LINE_COUNT lines)!"
            exit 1
        fi

        echo "✓ UniqueID added: \$LINE_COUNT lines (including header)"
        """
}
