process convertToFinalTSV {

    tag "$sample_id"

    input:
        tuple val(sample_id), path(annovar_txt)

    output:
        path("${sample_id}_Final_.tsv")

    script:
        """
        # Validate input file exists and is text
        if [ ! -f "${annovar_txt}" ]; then
            echo "ERROR: Input ANNOVAR file not found: ${annovar_txt}"
            exit 1
        fi

        if ! file "${annovar_txt}" | grep -qi "text"; then
            echo "ERROR: Input ANNOVAR file is corrupted (binary data detected)!"
            echo "File type: \$(file ${annovar_txt})"
            exit 1
        fi

        # Convert ANNOVAR multianno.txt to final TSV format with UniqueID
        # Input is already tab-separated from ANNOVAR

        awk 'BEGIN {FS=OFS="\\t"}
        NR==1 {
            # Header line: add UniqueID as first column
            print "UniqueID", \$0;
            next;
        }
        {
            # Data lines: create UniqueID as chr:pos:ref:alt (columns 1,2,4,5 in ANNOVAR output)
            uniqueID = \$1 ":" \$2 ":" \$4 ":" \$5;
            print uniqueID, \$0;
        }' ${annovar_txt} > ${sample_id}_Final_.tsv

        # Validate output
        if [ ! -s "${sample_id}_Final_.tsv" ]; then
            echo "ERROR: Output TSV file is empty!"
            exit 1
        fi

        # Check output is valid text
        if ! file "${sample_id}_Final_.tsv" | grep -qi "text"; then
            echo "ERROR: Output TSV file is corrupted!"
            exit 1
        fi

        # Verify we have at least a header line
        LINE_COUNT=\$(wc -l < "${sample_id}_Final_.tsv")
        if [ \$LINE_COUNT -lt 1 ]; then
            echo "ERROR: Output TSV has no content!"
            exit 1
        fi

        echo "TSV conversion completed: \$LINE_COUNT lines (including header)"
        """
}
