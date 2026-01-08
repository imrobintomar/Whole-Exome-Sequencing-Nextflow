process convertToFinalTSV {

    tag "$sample_id"

    input:
        tuple val(sample_id), path(annovar_txt)

    output:
        path("${sample_id}_Final_.tsv")

    script:
        """
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
        """
}
