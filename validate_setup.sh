#!/bin/bash

# =============================================================================
# WES Pipeline Pre-flight Validation Script
# =============================================================================
# This script validates your system setup before running the Nextflow pipeline
# Run this before executing the main pipeline to catch configuration issues
# =============================================================================

set -e  # Exit on error

ERRORS=0
WARNINGS=0

echo "╔══════════════════════════════════════════════════════════════╗"
echo "║   WES Pipeline Configuration Validation                     ║"
echo "╚══════════════════════════════════════════════════════════════╝"
echo ""

# =============================================================================
# 1. Check Required Software
# =============================================================================
echo "[1/7] Checking required software..."

check_software() {
    if command -v "$1" &> /dev/null; then
        VERSION=$($1 --version 2>&1 | head -n1 || echo "version unknown")
        echo "  ✓ $1 found: $VERSION"
    else
        echo "  ✗ $1 NOT FOUND"
        ((ERRORS++))
    fi
}

check_software bwa
check_software samtools
check_software gatk
check_software fastp
check_software bgzip
check_software tabix
check_software java

echo ""

# =============================================================================
# 2. Check Reference Files
# =============================================================================
echo "[2/7] Checking reference genome files..."

REF_DIR="/media/drprabudh/m1/hg38"
REF_FA="${REF_DIR}/hg38.fa"

if [[ -f "$REF_FA" ]]; then
    SIZE=$(du -h "$REF_FA" | cut -f1)
    echo "  ✓ Reference FASTA found: $REF_FA ($SIZE)"
else
    echo "  ✗ Reference FASTA NOT FOUND: $REF_FA"
    ((ERRORS++))
fi

if [[ -f "${REF_FA}.fai" ]]; then
    echo "  ✓ Reference index found: ${REF_FA}.fai"
else
    echo "  ✗ Reference index NOT FOUND: ${REF_FA}.fai"
    echo "    → Run: samtools faidx $REF_FA"
    ((ERRORS++))
fi

if [[ -f "${REF_DIR}/hg38.dict" ]]; then
    echo "  ✓ Reference dictionary found"
else
    echo "  ✗ Reference dictionary NOT FOUND: ${REF_DIR}/hg38.dict"
    echo "    → Run: gatk CreateSequenceDictionary -R $REF_FA"
    ((ERRORS++))
fi

echo ""

# =============================================================================
# 3. Check BWA Index
# =============================================================================
echo "[3/7] Checking BWA index files..."

BWA_INDEX_FILES=("amb" "ann" "bwt" "pac" "sa")
BWA_INDEX_PRESENT=true

for ext in "${BWA_INDEX_FILES[@]}"; do
    # Check for both hg38.ext and hg38.fa.ext patterns
    if [[ -f "${REF_DIR}/hg38.${ext}" ]] || [[ -f "${REF_DIR}/hg38.fa.${ext}" ]]; then
        if [[ -f "${REF_DIR}/hg38.${ext}" ]]; then
            echo "  ✓ hg38.${ext} found"
        else
            echo "  ✓ hg38.fa.${ext} found"
        fi
    else
        echo "  ✗ hg38.${ext} or hg38.fa.${ext} NOT FOUND"
        BWA_INDEX_PRESENT=false
        ((ERRORS++))
    fi
done

if [[ "$BWA_INDEX_PRESENT" == "false" ]]; then
    echo "    → Run: bwa index $REF_FA"
fi

echo ""

# =============================================================================
# 4. Check Known Sites VCF Files
# =============================================================================
echo "[4/7] Checking known sites VCF files..."

VCF_DIR="/media/drprabudh/m1/vcf_file"
KNOWN_SITES=(
    "Homo_sapiens_assembly38.known_indels.vcf.gz"
    "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
)

for vcf in "${KNOWN_SITES[@]}"; do
    VCF_PATH="${VCF_DIR}/${vcf}"
    if [[ -f "$VCF_PATH" ]]; then
        SIZE=$(du -h "$VCF_PATH" | cut -f1)
        echo "  ✓ $vcf ($SIZE)"

        # Check if indexed
        if [[ -f "${VCF_PATH}.tbi" ]]; then
            echo "    ✓ Index found"
        else
            echo "    ⚠ Index NOT FOUND: ${VCF_PATH}.tbi"
            echo "      → Run: tabix -p vcf $VCF_PATH"
            ((WARNINGS++))
        fi
    else
        echo "  ✗ $vcf NOT FOUND"
        ((ERRORS++))
    fi
done

echo ""

# =============================================================================
# 5. Check ANNOVAR Installation
# =============================================================================
echo "[5/7] Checking ANNOVAR installation..."

ANNOVAR_DIR="/media/drprabudh/m1/annovar"
ANNOVAR_DB="${ANNOVAR_DIR}/hg38_humandb"

if [[ -f "${ANNOVAR_DIR}/table_annovar.pl" ]]; then
    echo "  ✓ ANNOVAR script found"
else
    echo "  ✗ table_annovar.pl NOT FOUND in $ANNOVAR_DIR"
    ((ERRORS++))
fi

if [[ -d "$ANNOVAR_DB" ]]; then
    DB_COUNT=$(ls -1 "$ANNOVAR_DB" 2>/dev/null | wc -l)
    echo "  ✓ ANNOVAR database directory found ($DB_COUNT files)"

    # Check for required databases
    REQUIRED_DBS=("refGeneWithVer" "dbnsfp42a" "gnomad40_exome")
    for db in "${REQUIRED_DBS[@]}"; do
        if ls "$ANNOVAR_DB"/*"$db"* 1> /dev/null 2>&1; then
            echo "    ✓ $db database found"
        else
            echo "    ⚠ $db database NOT FOUND"
            ((WARNINGS++))
        fi
    done
else
    echo "  ✗ ANNOVAR database directory NOT FOUND: $ANNOVAR_DB"
    ((ERRORS++))
fi

echo ""

# =============================================================================
# 6. Check 1000 Genomes Files
# =============================================================================
echo "[6/7] Checking 1000 Genomes reference files..."

TG_DIR="/media/drprabudh/m1/annovar/1000"

if [[ -d "$TG_DIR" ]]; then
    echo "  ✓ 1000 Genomes directory found"

    # Check for chromosome files
    CHR_COUNT=0
    for chr in {1..22} X; do
        if [[ -f "${TG_DIR}/filtered_chr${chr}.vcf" ]] || [[ -f "${TG_DIR}/filtered_chr${chr}.vcf.gz" ]]; then
            ((CHR_COUNT++))
        fi
    done

    echo "    Found $CHR_COUNT/23 chromosome files"

    if [[ $CHR_COUNT -eq 23 ]]; then
        echo "    ✓ All chromosomes present"
    elif [[ $CHR_COUNT -gt 0 ]]; then
        echo "    ⚠ Some chromosomes missing"
        ((WARNINGS++))
    else
        echo "    ✗ No chromosome files found"
        ((ERRORS++))
    fi
else
    echo "  ✗ 1000 Genomes directory NOT FOUND: $TG_DIR"
    ((ERRORS++))
fi

echo ""

# =============================================================================
# 7. Check SnpSift JAR
# =============================================================================
echo "[7/7] Checking SnpSift installation..."

# Try to find SnpSift
SNPSIFT_LOCATIONS=(
    "/usr/local/share/snpsift-4.3.1t-2/SnpSift.jar"
    "/usr/share/java/snpsift.jar"
    "/opt/snpEff/SnpSift.jar"
)

SNPSIFT_FOUND=false
for jar in "${SNPSIFT_LOCATIONS[@]}"; do
    if [[ -f "$jar" ]]; then
        echo "  ✓ SnpSift JAR found: $jar"
        SNPSIFT_FOUND=true

        # Check if path matches config
        if grep -q "$jar" nextflow.config; then
            echo "    ✓ Path matches nextflow.config"
        else
            echo "    ⚠ Path differs from nextflow.config"
            echo "      → Update line 41 in nextflow.config to: $jar"
            ((WARNINGS++))
        fi
        break
    fi
done

if [[ "$SNPSIFT_FOUND" == "false" ]]; then
    echo "  ⚠ SnpSift JAR not found in common locations"
    echo "    Searching system..."
    FOUND_JAR=$(find /usr /opt -name "SnpSift.jar" 2>/dev/null | head -n1 || echo "")

    if [[ -n "$FOUND_JAR" ]]; then
        echo "    ✓ Found at: $FOUND_JAR"
        echo "      → Update line 41 in nextflow.config to: $FOUND_JAR"
        ((WARNINGS++))
    else
        echo "    ✗ SnpSift not found on system"
        echo "      → Install SnpEff/SnpSift package"
        ((ERRORS++))
    fi
fi

echo ""

# =============================================================================
# 8. Check Input FASTQ Files
# =============================================================================
echo "[Bonus] Checking input FASTQ files..."

INPUT_DIR="/media/drprabudh/m3/PRJNA855946/FASTQ"

if [[ -d "$INPUT_DIR" ]]; then
    FASTQ_COUNT=$(ls -1 "$INPUT_DIR"/*_{1,2}.fastq.gz 2>/dev/null | wc -l)
    PAIR_COUNT=$((FASTQ_COUNT / 2))

    if [[ $FASTQ_COUNT -gt 0 ]]; then
        echo "  ✓ Input directory found: $FASTQ_COUNT FASTQ files ($PAIR_COUNT pairs)"
    else
        echo "  ⚠ No FASTQ files found matching pattern *_{1,2}.fastq.gz"
        ((WARNINGS++))
    fi
else
    echo "  ⚠ Input directory not found: $INPUT_DIR"
    echo "    (This is OK if you'll specify a different path)"
    ((WARNINGS++))
fi

echo ""

# =============================================================================
# Summary
# =============================================================================
echo "════════════════════════════════════════════════════════════════"
echo "VALIDATION SUMMARY"
echo "════════════════════════════════════════════════════════════════"

if [[ $ERRORS -eq 0 ]] && [[ $WARNINGS -eq 0 ]]; then
    echo "✓ ALL CHECKS PASSED"
    echo ""
    echo "Your system is ready to run the WES pipeline!"
    echo ""
    echo "To start the pipeline:"
    echo "  nextflow run main.nf"
    exit 0
elif [[ $ERRORS -eq 0 ]]; then
    echo "⚠ VALIDATION PASSED WITH WARNINGS"
    echo ""
    echo "Errors:   $ERRORS"
    echo "Warnings: $WARNINGS"
    echo ""
    echo "The pipeline should run, but you may want to address warnings."
    exit 0
else
    echo "✗ VALIDATION FAILED"
    echo ""
    echo "Errors:   $ERRORS"
    echo "Warnings: $WARNINGS"
    echo ""
    echo "Please fix the errors above before running the pipeline."
    exit 1
fi
