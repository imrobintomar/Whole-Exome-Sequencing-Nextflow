# Pipeline Validation Summary

## ✅ Validation Results

Based on the validation run, here's the status of your system:

### 1. Software - ALL FOUND ✓
- ✅ bwa - found
- ✅ samtools 1.13 - found
- ✅ GATK 4.6.2.0 - found
- ✅ fastp 0.20.1 - found
- ✅ bgzip (htslib) 1.23 - found
- ✅ tabix (htslib) 1.23 - found
- ✅ java (openjdk 17.0.17) - found

### 2. Reference Genome - ALL PRESENT ✓
- ✅ hg38.fa (3.1G) - found
- ✅ hg38.fa.fai - found
- ✅ hg38.dict - found

### 3. BWA Index Files - ALL PRESENT ✓
- ✅ hg38.fa.amb - found
- ✅ hg38.fa.ann - found
- ✅ hg38.fa.bwt - found
- ✅ hg38.fa.pac - found
- ✅ hg38.fa.sa - found

### 4. Known Sites VCF Files - ALL PRESENT ✓
- ✅ Homo_sapiens_assembly38.known_indels.vcf.gz (59M) + index
- ✅ Mills_and_1000G_gold_standard.indels.hg38.vcf.gz (20M) + index

### 5. ANNOVAR Installation - COMPLETE ✓
- ✅ table_annovar.pl - found
- ✅ Database directory (90 files) - found
- ✅ refGeneWithVer database - found
- ✅ dbnsfp42a database - found
- ✅ gnomad40_exome database - found

### 6. 1000 Genomes Files - ALL PRESENT ✓
- ✅ Directory found: /media/drprabudh/m1/annovar/1000
- ✅ All 23 chromosomes present (chr1-22, chrX)

### 7. SnpSift JAR - NEEDS CONFIGURATION ⚠️
- ⚠️ Need to update nextflow.config with actual path

---

## 🔧 Required Action

### Update SnpSift Path in Configuration

The only remaining item is to configure the SnpSift JAR path.

**Find SnpSift on your system:**
```bash
# Search in conda environment (most likely location)
find ~/miniconda3 -name "SnpSift.jar" 2>/dev/null

# Or search common locations
find /usr/local /opt ~/tools -name "SnpSift.jar" 2>/dev/null
```

**Then update** [`nextflow.config`](nextflow.config) line 41:
```bash
nano nextflow.config
# Change line 41 to the actual path found above
```

---

## 📋 Input Data Check

### FASTQ Files

Check your input directory:
```bash
ls -lh /media/drprabudh/m3/PRJNA855946/FASTQ/*_{1,2}.fastq.gz
```

Expected format:
- SAMPLE1_1.fastq.gz
- SAMPLE1_2.fastq.gz
- SAMPLE2_1.fastq.gz
- SAMPLE2_2.fastq.gz

---

## ✨ System Status

### Overall: READY TO RUN (After SnpSift Configuration)

**Errors:** 0
**Warnings:** 1 (SnpSift path needs update)

All critical components are in place:
- ✅ All software installed
- ✅ Reference genome and indices complete
- ✅ Known sites VCF files indexed
- ✅ ANNOVAR fully configured
- ✅ 1000 Genomes data complete
- ⚠️ SnpSift path needs configuration

---

## 🚀 Next Steps

1. **Find and configure SnpSift:**
   ```bash
   find ~/miniconda3 -name "SnpSift.jar" 2>/dev/null
   nano nextflow.config  # Update line 41
   ```

2. **Verify your FASTQ files:**
   ```bash
   ls /media/drprabudh/m3/PRJNA855946/FASTQ/*_{1,2}.fastq.gz
   ```

3. **Run a test:**
   ```bash
   # Test with help command
   nextflow run main.nf --help
   ```

4. **Launch the pipeline:**
   ```bash
   # Make sure bgzip/tabix are in PATH
   export PATH="/home/drprabudh/miniconda3/envs/lncrna/bin:$PATH"

   # Run with full reporting
   nextflow run main.nf -with-trace -with-report -with-timeline -with-dag
   ```

---

## 💡 Tips

### Keep bgzip/tabix Available

For every session, either:

**Option A:** Activate conda environment
```bash
conda activate lncrna
```

**Option B:** Add to PATH
```bash
export PATH="/home/drprabudh/miniconda3/envs/lncrna/bin:$PATH"
```

**Option C:** Make permanent (add to ~/.bashrc)
```bash
echo 'export PATH="/home/drprabudh/miniconda3/envs/lncrna/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

### Disk Space

Ensure you have adequate space:
```bash
df -h /media/drprabudh/m3
# Need at least 500 GB free
```

### Expected Runtime

- Single sample: 24-48 hours
- Multiple samples: Processed in parallel

---

## 📞 If You Encounter Issues

1. Check PATH has bgzip/tabix: `which bgzip tabix`
2. Verify SnpSift path in config: `grep snpsift_jar nextflow.config`
3. Check FASTQ file naming: Must be `*_{1,2}.fastq.gz`
4. Review logs: `tail -f .nextflow.log`

---

**Your system is ready! Just configure SnpSift and you're good to go.**
