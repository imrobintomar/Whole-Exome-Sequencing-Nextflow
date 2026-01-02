#!/bin/bash

echo "╔════════════════════════════════════════════════════════════════╗"
echo "║                                                                ║"
echo "║  STARTING PIPELINE WITH STACKOVERFLOW FIX                      ║"
echo "║                                                                ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""
echo "Setting JVM options to prevent StackOverflowError..."
export NXF_OPTS="-Xms2g -Xmx8g -Xss4m"
echo "✓ NXF_OPTS = $NXF_OPTS"
echo ""
echo "Running pipeline with -resume flag..."
echo ""

nextflow run main.nf -resume "$@"

echo ""
echo "Pipeline execution completed."
