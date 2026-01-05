#!/bin/bash

# Download gnomAD gene constraint metrics for ACMG classification
# This provides pLI (probability of LOF intolerance) and LOEUF scores

echo "📥 Downloading gnomAD v4 gene constraint metrics..."

# Create data directory
mkdir -p data

# Download gnomAD constraint file (~3MB)
curl -L "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv" \
  -o data/gnomad_constraint.tsv

if [ $? -eq 0 ]; then
    echo "✅ Downloaded successfully to data/gnomad_constraint.tsv"
    echo "📊 File size: $(du -h data/gnomad_constraint.tsv | cut -f1)"
    echo "📈 Gene count: $(wc -l < data/gnomad_constraint.tsv)"
else
    echo "❌ Download failed"
    exit 1
fi

echo ""
echo "Sample data:"
head -n 3 data/gnomad_constraint.tsv

echo ""
echo "✅ Setup complete! You can now use ACMG classification."
