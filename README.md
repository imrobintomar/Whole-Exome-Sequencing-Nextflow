# Whole Exome Sequencing (WES) Pipeline

Production-ready Nextflow DSL2 pipeline for processing whole exome sequencing data from raw FASTQ files to annotated, filtered variants.

## Pipeline Overview

```
FASTQ Files
    ↓
[1] QC & Filtering (fastp)
    ↓
[2] Alignment (BWA-MEM)
    ↓
[3] Sorting & QC (GATK SortSam + samtools flagstat)
    ↓
[4] Mark Duplicates (GATK MarkDuplicates)
    ↓
[5] Base Quality Score Recalibration (GATK BQSR)
    ↓
[6] Variant Calling (GATK HaplotypeCaller)
    ↓
[7] Annotation (1000 Genomes + ANNOVAR)
    ↓
[8] Filtering (SnpSift)
    ↓
Final TSV Output
```

## Recent Fixes (v1.1)

All critical issues have been resolved:
- ✓ Fixed process name mismatches
- ✓ Added BAM indexing for GATK compatibility
- ✓ Removed hardcoded resources for dynamic allocation
- ✓ Fixed filtering column numbers (AF and DP)
- ✓ Added VCF compression and indexing
- ✓ Added error handling for missing 1000G files
- ✓ Replaced hardcoded thread counts with task.cpus

See [CHANGELOG.md](CHANGELOG.md) for detailed changes.

## Prerequisites

### Software Requirements

| Tool | Minimum Version | Purpose |
|------|----------------|---------|
| Nextflow | ≥21.04 | Workflow engine |
| BWA | ≥0.7.17 | Read alignment |
| SAMtools | ≥1.17 | BAM manipulation |
| GATK | ≥4.6.0 | Variant calling & processing |
| fastp | ≥0.23.0 | Quality control |
| bgzip/tabix | (any) | VCF compression/indexing |
| Java | ≥1.8 | Required for GATK/SnpSift |
| ANNOVAR | (latest) | Variant annotation |
| SnpEff/SnpSift | ≥4.3 | Variant filtering |

### Reference Data Requirements

1. **Human Reference Genome (hg38)**
   - `hg38.fa` - FASTA file
   - `hg38.fa.fai` - FASTA index (create with `samtools faidx`)
   - `hg38.dict` - Sequence dictionary (create with `gatk CreateSequenceDictionary`)
   - BWA index files (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`)

2. **Known Sites VCF Files**
   - Homo_sapiens_assembly38.known_indels.vcf.gz
   - Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
   - Both must be indexed with tabix (`.tbi` files)

3. **ANNOVAR Databases**
   - refGeneWithVer
   - dbnsfp42a
   - clinvar_20240416
   - gnomad40_exome
   - avsnp150
   - cosmic84_coding
   - exac03

4. **1000 Genomes VCF Files**
   - `filtered_chr{1..22,X}.vcf` or `.vcf.gz`

## Installation

### 1. Clone/Download Pipeline

```bash
cd /path/to/your/workspace
# Pipeline files should be in ./WholeExome/
```

### 2. Install Dependencies

**Option A: Using Conda (Recommended)**
```bash
conda create -n wes-pipeline -c bioconda \
    nextflow bwa samtools gatk4 fastp bcftools
conda activate wes-pipeline
```

**Option B: System Package Manager**
```bash
# Ubuntu/Debian
sudo apt-get install bwa samtools bcftools

# Download GATK from https://github.com/broadinstitute/gatk/releases
# Install fastp from https://github.com/OpenGene/fastp
```

### 3. Configure Paths

Edit [`nextflow.config`](nextflow.config) to match your system:

```groovy
params {
    // INPUT/OUTPUT
    input_dir   = '/path/to/your/fastq/files'
    output_dir  = '/path/to/output/directory'

    // REFERENCE GENOME
    reference        = '/path/to/hg38/hg38.fa'
    reference_index  = '/path/to/hg38/hg38.fa.fai'
    reference_dict   = '/path/to/hg38/hg38.dict'
    bwa_index        = '/path/to/hg38/hg38'

    // KNOWN SITES
    known_sites = [
        '/path/to/Homo_sapiens_assembly38.known_indels.vcf.gz',
        '/path/to/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
    ]

    // ANNOTATION TOOLS
    annovar_dir      = '/path/to/annovar'
    annovar_db       = '/path/to/annovar/hg38_humandb'
    annovar_xreffile = '/path/to/annovar/dbnsfp/gene/dbNSFP4.7_gene'
    thousand_genomes_dir = '/path/to/1000genomes/vcfs'
    snpsift_jar      = '/path/to/SnpSift.jar'

    // FILTERING THRESHOLDS
    max_af     = 0.05   // Maximum allele frequency
    min_depth = 5       // Minimum read depth
}
```

### 4. Validate Setup

Run the validation script to check your configuration:

```bash
cd /media/drprabudh/m3/Nextflow-Script/WholeExome
./validate_setup.sh
```

This will check:
- Software availability
- Reference file existence
- Index file presence
- Database completeness
- Configuration correctness

## Usage

### Basic Execution

```bash
nextflow run main.nf
```

### Common Options

```bash
# Custom input directory
nextflow run main.nf --input_dir /path/to/fastq

# Custom output directory
nextflow run main.nf --output_dir /path/to/results

# Resume failed pipeline
nextflow run main.nf -resume

# Generate execution reports
nextflow run main.nf -with-trace -with-report -with-timeline -with-dag
```

### Resource Management Profiles

Control resource allocation:

```bash
# Default: 80% of system resources
nextflow run main.nf -profile standard

# Conservative: 50% of resources (for shared systems)
nextflow run main.nf -profile conservative

# Aggressive: 95% of resources (dedicated systems)
nextflow run main.nf -profile aggressive
```

### Container Execution (Future)

```bash
# Using Docker
nextflow run main.nf -profile docker

# Using Singularity
nextflow run main.nf -profile singularity
```

### Help Message

```bash
nextflow run main.nf --help
```

## Input Data Format

The pipeline expects paired-end FASTQ files with this naming convention:

```
SAMPLE1_1.fastq.gz
SAMPLE1_2.fastq.gz
SAMPLE2_1.fastq.gz
SAMPLE2_2.fastq.gz
```

Pattern: `*_{1,2}.fastq.gz`

Place all FASTQ files in the input directory specified by `params.input_dir`.

## Output Structure

```
output_dir/
├── filtered_fastp/          # QC-filtered FASTQ files
│   ├── SAMPLE1_1_filtered.fastq.gz
│   ├── SAMPLE1_2_filtered.fastq.gz
│   └── SAMPLE1.html         # QC report
├── Mapsam/                  # BAM files
│   ├── SAMPLE1.sorted.bam
│   ├── SAMPLE1.sorted.bam.bai
│   ├── SAMPLE1_markdup.sorted.bam
│   ├── SAMPLE1_markdup.sorted.bam.bai
│   ├── SAMPLE1_recall.bam   # Final BAM
│   ├── SAMPLE1.Stat.txt     # Flagstat report
│   └── SAMPLE1_recall.table # BQSR recalibration table
├── Germline_VCF/            # Variant files
│   ├── SAMPLE1.vcf.gz
│   ├── SAMPLE1.vcf.gz.tbi
│   ├── SAMPLE1_1000genome.vcf.gz
│   ├── SAMPLE1.annovar.hg38_multianno.vcf
│   └── SAMPLE1_final_annotated.tsv  # FINAL OUTPUT
└── logs/                    # Execution reports
    ├── trace.txt
    ├── report.html
    ├── timeline.html
    └── pipeline.svg         # DAG visualization
```

## Output Files Explained

### Final Output: `*_final_annotated.tsv`

Tab-separated file containing filtered variants with annotations:

| Column | Description |
|--------|-------------|
| CHROM | Chromosome |
| POS | Position |
| ID | Variant ID (dbSNP) |
| REF | Reference allele |
| ALT | Alternate allele |
| QUAL | Variant quality score |
| FILTER | Filter status |
| DP | Read depth |
| AC | Allele count |
| AF | Allele frequency |
| AN | Total alleles |
| ExonicFunc | Functional effect |
| Gene | Gene name |
| SIFT_pred | SIFT prediction |
| Polyphen2_HDIV_pred | PolyPhen-2 prediction |
| ... | Additional annotations |

**Filtering Applied:**
- Allele frequency (AF) ≤ 0.05
- Read depth (DP) ≥ 5
- Non-missing values

## Troubleshooting

### Pipeline Fails at Alignment

**Symptom:** BWA-MEM process fails
**Solution:** Check BWA index files exist and are complete
```bash
ls -lh /path/to/hg38/hg38.{amb,ann,bwt,pac,sa}
```

### GATK Tools Fail with "Index not found"

**Symptom:** BaseRecalibrator or HaplotypeCaller fails
**Solution:** The updated pipeline automatically creates BAM indices. If error persists:
```bash
# Manually index a BAM file
samtools index your_file.bam
```

### Known Sites Files Not Found

**Symptom:** BaseRecalibrator can't find known sites
**Solution:** Check VCF files are indexed
```bash
tabix -p vcf /path/to/known_sites.vcf.gz
```

### Out of Memory Errors

**Symptom:** Process killed with exit code 137
**Solution:** The pipeline has automatic retry with increased memory. If still failing:
1. Use conservative profile: `-profile conservative`
2. Increase `maxMemoryGB` in [nextflow.config:52](nextflow.config#L52)

### SnpSift Not Found

**Symptom:** Annotation or filtering fails
**Solution:** Update `params.snpsift_jar` in config:
```bash
# Find SnpSift
find /usr /opt -name "SnpSift.jar"

# Update nextflow.config line 41
params.snpsift_jar = '/actual/path/to/SnpSift.jar'
```

### 1000 Genomes Annotation Warnings

**Symptom:** "Warning: 1000G file for chrX not found"
**Solution:** This is informational - pipeline continues with available chromosomes. To fix:
- Download missing chromosome VCFs
- Ensure naming: `filtered_chr{1..22,X}.vcf[.gz]`

## Performance Optimization

### Recommended Resources

| Step | CPUs | Memory | Time (per sample) |
|------|------|--------|-------------------|
| fastp | 4-8 | 8 GB | 30-60 min |
| BWA-MEM | 16-32 | 32-64 GB | 4-8 hours |
| GATK Sort | 4-8 | 16-32 GB | 1-2 hours |
| MarkDuplicates | 8-16 | 32-48 GB | 2-4 hours |
| BQSR | 16-32 | 32-48 GB | 4-8 hours |
| HaplotypeCaller | 16-32 | 32-64 GB | 6-12 hours |
| Annotation | 8-16 | 16-32 GB | 2-6 hours |

**Total Runtime:** ~24-48 hours per sample on a typical workstation

### Speeding Up the Pipeline

1. **Increase parallelism:** Process multiple samples simultaneously
2. **Use faster storage:** SSD/NVMe for working directory
3. **Optimize resources:** Use `-profile aggressive` on dedicated systems
4. **Enable caching:** Always use `-resume` when rerunning

## Advanced Configuration

### Custom Filtering Thresholds

Edit [nextflow.config:43-44](nextflow.config#L43-L44):
```groovy
params.max_af     = 0.01   // Stricter: rare variants only
params.min_depth = 10      // Higher quality threshold
```

### Modify Resource Allocation

Edit resource percentages in [nextflow.config:54-58](nextflow.config#L54-L58):
```groovy
def modeScale = [
    auto:         [cpu: 0.8, mem: 0.8],   // Default
    conservative: [cpu: 0.5, mem: 0.5],   // Safer
    aggressive:   [cpu: 0.95, mem: 0.95]  // Maximum
]
```

### Process-Specific Settings

Adjust individual process resources in [nextflow.config:96-187](nextflow.config#L96-L187):
```groovy
withName: bwaMem {
    cpus   = cpusPct(90)   // Use 90% of available CPUs
    memory = memPct(80)    // Use 80% of available memory
    time   = '24h'         // Max runtime
}
```

## Citation

If you use this pipeline, please cite the tools:

- **Nextflow:** Di Tommaso, P., et al. (2017). Nat Biotechnol. 35, 316-319
- **BWA:** Li, H. and Durbin, R. (2009). Bioinformatics 25, 1754-1760
- **GATK:** McKenna, A., et al. (2010). Genome Res. 20, 1297-1303
- **fastp:** Chen, S., et al. (2018). Bioinformatics 34, i884-i890
- **ANNOVAR:** Wang, K., et al. (2010). Nucleic Acids Res. 38, e164

## Support

- **Issues:** Report bugs or request features via GitHub issues
- **Documentation:** See [CHANGELOG.md](CHANGELOG.md) for version history
- **Validation:** Run `./validate_setup.sh` before reporting issues

## License

This pipeline is provided as-is for academic and research use.

## Author

Robin Tomar
Version: 1.1
Last Updated: 2026-01-01
