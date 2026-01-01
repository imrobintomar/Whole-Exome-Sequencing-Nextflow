# Pre-Run Checklist

Complete this checklist before running the pipeline for the first time.

## ☐ Phase 1: System Requirements

### Software Installation
- [ ] Nextflow installed and in PATH (`nextflow -version`)
- [ ] BWA installed (`bwa 2>&1 | head -3`)
- [ ] SAMtools installed (`samtools --version`)
- [ ] GATK installed (`gatk --version`)
- [ ] fastp installed (`fastp --version`)
- [ ] bgzip/tabix installed (`bgzip --version && tabix --version`)
- [ ] Java installed (`java -version`)

### Reference Data Files
- [ ] hg38.fa exists and is accessible
- [ ] hg38.fa.fai exists (or create: `samtools faidx hg38.fa`)
- [ ] hg38.dict exists (or create: `gatk CreateSequenceDictionary -R hg38.fa`)
- [ ] BWA index files exist (*.amb, *.ann, *.bwt, *.pac, *.sa)
  - If missing: `bwa index hg38.fa`

### Known Sites VCF Files
- [ ] Homo_sapiens_assembly38.known_indels.vcf.gz exists
- [ ] Homo_sapiens_assembly38.known_indels.vcf.gz.tbi exists
  - If missing: `tabix -p vcf Homo_sapiens_assembly38.known_indels.vcf.gz`
- [ ] Mills_and_1000G_gold_standard.indels.hg38.vcf.gz exists
- [ ] Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi exists
  - If missing: `tabix -p vcf Mills_and_1000G_gold_standard.indels.hg38.vcf.gz`

---

## ☐ Phase 2: Configuration

### Find SnpSift Installation
```bash
find /usr /opt -name "SnpSift.jar" 2>/dev/null
```
- [ ] Found SnpSift.jar at path: _______________________________
- [ ] Updated nextflow.config line 41 with correct path

### Verify ANNOVAR Setup
- [ ] table_annovar.pl exists at configured path
- [ ] ANNOVAR database directory exists (hg38_humandb)
- [ ] Required databases present:
  - [ ] refGeneWithVer
  - [ ] dbnsfp42a
  - [ ] gnomad40_exome
  - [ ] clinvar_20240416

### Check 1000 Genomes Files
- [ ] 1000 Genomes directory exists
- [ ] Chromosome files present (filtered_chr{1..22,X}.vcf or .vcf.gz)
- [ ] At least 20/23 chromosomes available

---

## ☐ Phase 3: Input Data

### FASTQ Files
```bash
ls -lh /media/drprabudh/m3/PRJNA855946/FASTQ/*_{1,2}.fastq.gz
```
- [ ] FASTQ files exist in input directory
- [ ] Files follow naming pattern: SAMPLENAME_{1,2}.fastq.gz
- [ ] Each sample has both _1 and _2 files (paired-end)
- [ ] Number of sample pairs: _______

### Disk Space
```bash
df -h /media/drprabudh/m3
```
- [ ] At least 500 GB free space available
- [ ] Current free space: _______ GB

---

## ☐ Phase 4: Validation

### Run Setup Validation
```bash
cd /media/drprabudh/m3/Nextflow-Script/WholeExome
./validate_setup.sh
```
- [ ] Validation script runs successfully
- [ ] Zero ERRORS reported
- [ ] All WARNINGS reviewed and understood

### Test Pipeline Syntax
```bash
nextflow run main.nf --help
```
- [ ] Help message displays correctly
- [ ] No syntax errors

---

## ☐ Phase 5: Configuration Verification

### Verify Paths in nextflow.config
Open `nextflow.config` and verify these paths:

**Line 24:** `params.input_dir`
- [ ] Path exists: _______________________________
- [ ] Contains FASTQ files

**Line 25:** `params.output_dir`
- [ ] Parent directory exists
- [ ] Has write permissions

**Line 27:** `params.reference`
- [ ] File exists: /media/drprabudh/m1/hg38/hg38.fa

**Line 28:** `params.reference_index`
- [ ] File exists: /media/drprabudh/m1/hg38/hg38.fa.fai

**Line 29:** `params.reference_dict`
- [ ] File exists: /media/drprabudh/m1/hg38/hg38.dict

**Line 30:** `params.bwa_index`
- [ ] All BWA index files exist with this prefix

**Lines 32-35:** `params.known_sites`
- [ ] Both VCF files exist
- [ ] Both have .tbi index files

**Line 37:** `params.annovar_dir`
- [ ] Directory exists
- [ ] Contains table_annovar.pl

**Line 38:** `params.annovar_db`
- [ ] Directory exists
- [ ] Contains database files

**Line 40:** `params.thousand_genomes_dir`
- [ ] Directory exists
- [ ] Contains filtered_chr*.vcf files

**Line 41:** `params.snpsift_jar`
- [ ] File exists
- [ ] Correct path from Phase 2

---

## ☐ Phase 6: Resource Planning

### System Resources
```bash
nproc  # Number of CPU cores
free -h  # Available memory
```
- [ ] CPU cores available: _______
- [ ] Total RAM: _______ GB
- [ ] Free RAM: _______ GB

### Choose Resource Profile
- [ ] **Standard** (default, 80% resources) - Recommended
- [ ] **Conservative** (50% resources) - For shared systems
- [ ] **Aggressive** (95% resources) - For dedicated systems

Selected profile: _______________________________

---

## ☐ Phase 7: Test Run (OPTIONAL BUT RECOMMENDED)

### Single Sample Test
```bash
# Backup other samples
mkdir -p ../fastq_backup
cp /media/drprabudh/m3/PRJNA855946/FASTQ/* ../fastq_backup/

# Keep only one sample pair in input directory
# Delete or move others temporarily

# Run pipeline on single sample
nextflow run main.nf -with-trace -with-report

# Check if successful
ls -lh /media/drprabudh/m3/PRJNA855946/Germline_VCF/*_final_annotated.tsv

# Restore all samples if successful
cp ../fastq_backup/* /media/drprabudh/m3/PRJNA855946/FASTQ/
```

Test run results:
- [ ] Pipeline completed without errors
- [ ] Final TSV file generated
- [ ] File contains variants (not just header)

---

## ☐ Phase 8: Ready to Launch

### Final Checks
- [ ] All previous phases completed
- [ ] No critical errors from validation
- [ ] Sufficient disk space
- [ ] Test run successful (if performed)
- [ ] Have reviewed expected runtime (~24-48 hours per sample)

### Backup Strategy
- [ ] Important data backed up (input FASTQ files)
- [ ] Understand that work/ directory will use ~200-300 GB
- [ ] Plan to run `nextflow clean` after successful completion

### Monitoring Plan
- [ ] Know how to monitor: `tail -f .nextflow.log`
- [ ] Know where reports will be: `/media/drprabudh/m3/PRJNA855946/logs/`
- [ ] Understand resume functionality: `nextflow run main.nf -resume`

---

## ☐ Launch Commands

### First Full Run
```bash
cd /media/drprabudh/m3/Nextflow-Script/WholeExome

# With full reporting
nextflow run main.nf -with-trace -with-report -with-timeline -with-dag

# OR with custom profile
nextflow run main.nf -profile conservative -with-trace -with-report
```

### If Pipeline Fails/Stops
```bash
# Always use -resume to avoid restarting from scratch
nextflow run main.nf -resume -with-trace -with-report
```

---

## Common Issues Reference

| Issue | Check | Solution |
|-------|-------|----------|
| "No FASTQ files found" | FASTQ naming | Ensure *_{1,2}.fastq.gz pattern |
| "Reference not found" | Path in config | Update params.reference in config |
| Exit code 137 | Memory | Use `-profile conservative` |
| "Index not found" | .bai files | Fixed in new version (auto-creates) |
| "SnpSift not found" | JAR path | Update line 41 in nextflow.config |

---

## Post-Run Checklist

After pipeline completes successfully:

- [ ] Verify all samples have output TSV files
- [ ] Check variant counts are reasonable (not 0, not millions)
- [ ] Review QC reports (fastp HTML files)
- [ ] Check alignment statistics (.Stat.txt files)
- [ ] Review execution reports (report.html, timeline.html)
- [ ] Backup important results
- [ ] Clean work directory: `rm -rf work/` (saves ~200-300 GB)

---

## Sign-Off

I certify that I have completed all applicable items in this checklist:

- **Date:** _______________________________
- **Completed by:** _______________________________
- **Notes:** _______________________________
         _______________________________
         _______________________________

---

## Quick Reference

**Validate:** `./validate_setup.sh`
**Help:** `nextflow run main.nf --help`
**Run:** `nextflow run main.nf -with-trace -with-report`
**Resume:** `nextflow run main.nf -resume`
**Monitor:** `tail -f .nextflow.log`
**Reports:** `firefox /media/drprabudh/m3/PRJNA855946/logs/report.html`

**Documentation:**
- [README.md](README.md) - Full documentation
- [QUICK_START.md](QUICK_START.md) - Quick reference
- [FIXES_SUMMARY.md](FIXES_SUMMARY.md) - What was fixed
- [CHANGELOG.md](CHANGELOG.md) - Detailed changes
