process annovarToTSV {

    tag "$sample_id"

    publishDir "${params.outdir}/tsv", mode: 'copy'

    input:
        tuple val(sample_id), path(multianno_txt)

    output:
        tuple val(sample_id), path("${sample_id}.annovar.tsv")

    errorStrategy 'terminate'
    maxRetries 1

    script:
    """
    echo "Converting ANNOVAR → TSV for ${sample_id}"
    echo "Input: ${multianno_txt}"

    convert_annovar_to_tsv.py \
        ${multianno_txt} \
        ${sample_id}.annovar.tsv

    if [ ! -s ${sample_id}.annovar.tsv ]; then
        echo "ERROR: TSV output is empty!"
        exit 1
    fi

    echo "TSV conversion completed"
    """
}

process addUniqueID {

    tag "$sample_id"

    publishDir "${params.output_dir}/Final_Annotated", mode: 'copy'

    input:
        tuple val(sample_id), path(annovar_txt)

    output:
        tuple val(sample_id), path("${sample_id}_Final_.txt")

    errorStrategy 'terminate'
    maxRetries 1

    script:
    def script_path = "${projectDir}/processes/convert_annovar_to_tsv.py"
    """
    # Validate input file exists
    if [ ! -f "${annovar_txt}" ]; then
        echo "ERROR: Input ANNOVAR file not found: ${annovar_txt}"
        exit 1
    fi

    # Add UniqueID column to ANNOVAR output
    python3 ${script_path} \
        "${annovar_txt}" \
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
