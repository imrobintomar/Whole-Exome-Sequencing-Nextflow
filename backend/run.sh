#!/bin/bash

# ==========================================
# Quick Pipeline Runner with JVM Fix
# ==========================================
# This fixes the StackOverflowError by setting proper JVM memory
# Run this instead of: nextflow run ../main.nf
# ==========================================

# Set JVM options to fix StackOverflowError
export NXF_OPTS="-Xms2g -Xmx8g -Xss4m"

# Change to parent directory where main.nf is located
cd "$(dirname "$0")/.."

# Run pipeline with all arguments passed through
nextflow run main.nf "$@"
