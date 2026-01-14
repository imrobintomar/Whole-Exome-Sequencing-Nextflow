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
