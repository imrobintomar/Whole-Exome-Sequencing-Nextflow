# WES Pipeline - Quick Start Guide

## Pre-Flight Checklist

```bash
# 1. Navigate to pipeline directory
cd /media/drprabudh/m3/Nextflow-Script/WholeExome

# 2. Validate your setup
./validate_setup.sh

# 3. Check your FASTQ files are named correctly
ls -lh /media/drprabudh/m3/PRJNA855946/FASTQ/*_{1,2}.fastq.gz
```

## First Time Setup

### 1. Find and Configure SnpSift Path

```bash
# Find SnpSift.jar on your system
find /usr /opt -name "SnpSift.jar" 2>/dev/null

# Update nextflow.config line 41 with the actual path
nano nextflow.config
# Change: params.snpsift_jar = '/actual/path/to/SnpSift.jar'
```

### 2. Create Missing Indices (if needed)

```bash
# Reference FASTA index
samtools faidx /media/drprabudh/m1/hg38/hg38.fa

# Reference dictionary
gatk CreateSequenceDictionary -R /media/drprabudh/m1/hg38/hg38.fa

# BWA index
bwa index /media/drprabudh/m1/hg38/hg38.fa

# Index known sites VCF files
tabix -p vcf /media/drprabudh/m1/vcf_file/Homo_sapiens_assembly38.known_indels.vcf.gz
tabix -p vcf /media/drprabudh/m1/vcf_file/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
```

## Running the Pipeline

### Standard Run

```bash
# Default settings (uses 80% of system resources)
nextflow run main.nf
```

### With Full Reporting

```bash
# Recommended for first run
nextflow run main.nf -with-trace -with-report -with-timeline -with-dag
```

### Resume After Failure

```bash
# ALWAYS use -resume to avoid re-running completed steps
nextflow run main.nf -resume
```

### Custom Paths

```bash
# Override input/output directories
nextflow run main.nf \
    --input_dir /path/to/my/fastq \
    --output_dir /path/to/my/results
```

### Resource Profiles

```bash
# Shared system (uses 50% of resources)
nextflow run main.nf -profile conservative

# Dedicated system (uses 95% of resources)
nextflow run main.nf -profile aggressive
```

## Monitoring Progress

### Real-time Monitoring

```bash
# In another terminal, monitor the trace file
tail -f /media/drprabudh/m3/PRJNA855946/logs/trace.txt

# Or use watch to refresh every 5 seconds
watch -n 5 'tail -20 /media/drprabudh/m3/PRJNA855946/logs/trace.txt'
```

### Check Process Status

```bash
# View Nextflow log
tail -f .nextflow.log

# List work directories
ls -lhtr work/*/*
```

### View Reports (after completion)

```bash
# Open HTML report in browser
firefox /media/drprabudh/m3/PRJNA855946/logs/report.html

# View timeline
firefox /media/drprabudh/m3/PRJNA855946/logs/timeline.html

# View pipeline DAG
firefox /media/drprabudh/m3/PRJNA855946/logs/pipeline.svg
```

## Common Issues & Quick Fixes

### Issue: "No FASTQ files found"

```bash
# Check files exist and naming is correct
ls -lh /media/drprabudh/m3/PRJNA855946/FASTQ/*_{1,2}.fastq.gz

# Files should be named like: SAMPLE1_1.fastq.gz, SAMPLE1_2.fastq.gz
```

### Issue: "Reference not found"

```bash
# Verify reference exists
ls -lh /media/drprabudh/m1/hg38/hg38.fa

# Check config points to correct path
grep "params.reference" nextflow.config
```

### Issue: Process fails with exit code 137 (Out of Memory)

```bash
# Use conservative profile
nextflow run main.nf -profile conservative -resume

# Or edit nextflow.config line 52 to reduce maxMemoryGB
```

### Issue: "SnpSift.jar not found"

```bash
# Find SnpSift
find /usr /opt /home -name "SnpSift.jar" 2>/dev/null

# Update config with correct path
nano nextflow.config  # Line 41
```

### Issue: "known sites not found"

```bash
# Check VCF files exist
ls -lh /media/drprabudh/m1/vcf_file/*.vcf.gz

# Check they're indexed
ls -lh /media/drprabudh/m1/vcf_file/*.vcf.gz.tbi

# If missing index:
tabix -p vcf /media/drprabudh/m1/vcf_file/YOUR_FILE.vcf.gz
```

## Pipeline Cleanup

### Remove Work Directory (saves space)

```bash
# ONLY after successful completion
rm -rf work/

# Be careful! This removes all intermediate files
# You'll need to rerun from scratch if you delete this
```

### Clean Everything (start fresh)

```bash
# Remove work directory and Nextflow metadata
nextflow clean -f

# Remove all generated files
rm -rf work/ .nextflow/ .nextflow.log*
```

## Performance Tips

### Speed Up Multiple Runs

```bash
# Always use -resume to skip completed tasks
nextflow run main.nf -resume

# Nextflow will automatically detect what needs to be rerun
```

### Parallel Processing

```bash
# The pipeline automatically processes multiple samples in parallel
# No special flags needed - just put all FASTQ files in input_dir
```

### Use Fast Storage

```bash
# Set work directory to SSD/NVMe
export NXF_WORK="/path/to/fast/storage/work"
nextflow run main.nf
```

## What to Expect

### Timeline (per sample)

1. **fastp QC**: 30-60 minutes
2. **BWA alignment**: 4-8 hours
3. **Sorting**: 1-2 hours
4. **Mark duplicates**: 2-4 hours
5. **BQSR**: 4-8 hours
6. **Variant calling**: 6-12 hours
7. **Annotation**: 2-6 hours
8. **Filtering**: 5-10 minutes

**Total**: ~24-48 hours per sample

### Disk Space Requirements

- Input FASTQ (per sample): ~10-20 GB
- Work directory: ~200-300 GB (temporary)
- Final output: ~50-100 GB

**Recommendation**: Have at least 500 GB free space

### Final Output Location

```bash
# Main result file
/media/drprabudh/m3/PRJNA855946/Germline_VCF/SAMPLE_final_annotated.tsv

# This is your filtered, annotated variant list
```

## Getting Help

```bash
# Show pipeline help
nextflow run main.nf --help

# Validate setup
./validate_setup.sh

# Check Nextflow version
nextflow -version

# Check software versions
bwa 2>&1 | head -3
samtools --version
gatk --version
```

## After Pipeline Completes

### 1. Verify Output

```bash
# Check final TSV file exists and has content
wc -l /media/drprabudh/m3/PRJNA855946/Germline_VCF/*_final_annotated.tsv

# Should have variants (more than just header line)
```

### 2. Review QC Metrics

```bash
# Check alignment statistics
cat /media/drprabudh/m3/PRJNA855946/Mapsam/*.Stat.txt

# Review fastp HTML reports
firefox /media/drprabudh/m3/PRJNA855946/filtered_fastp/*.html
```

### 3. Analyze Results

```bash
# Count variants
tail -n +2 /media/drprabudh/m3/PRJNA855946/Germline_VCF/SAMPLE_final_annotated.tsv | wc -l

# View top variants
head -20 /media/drprabudh/m3/PRJNA855946/Germline_VCF/SAMPLE_final_annotated.tsv | column -t
```

## Next Steps

- Import TSV into Excel/R for downstream analysis
- Filter further based on pathogenicity predictions
- Prioritize variants based on clinical significance
- Validate candidates with Sanger sequencing

---

**Need More Help?**
- Check [README.md](README.md) for detailed documentation
- See [CHANGELOG.md](CHANGELOG.md) for what was fixed
- Run `./validate_setup.sh` to diagnose configuration issues
