#!/usr/bin/env python3
"""
Diagnostic script to check why jobs are failing
"""

import sqlite3
import os
from pathlib import Path

def diagnose_jobs():
    print("=" * 60)
    print("WES Pipeline Job Diagnostics")
    print("=" * 60)
    print()

    # Check database
    db_path = "wes_pipeline.db"
    if not os.path.exists(db_path):
        print(" Database file not found!")
        print("   The backend may not have been started yet.")
        return

    db_size = os.path.getsize(db_path)
    print(f" Database: {db_path} ({db_size} bytes)")

    if db_size == 0:
        print("  Database is empty (0 bytes)")
        print("   This means the backend hasn't initialized the database yet.")
        print()
        print(" Solution:")
        print("   1. Make sure the backend is running: python3 main.py")
        print("   2. The database will be initialized on first startup")
        return

    # Connect to database
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    try:
        # Get recent jobs
        cursor.execute('''
            SELECT job_id, sample_name, status, current_step,
                   error_message, created_at, fastq_r1_path, fastq_r2_path
            FROM jobs
            ORDER BY created_at DESC
            LIMIT 5
        ''')
        jobs = cursor.fetchall()

        if not jobs:
            print(" No jobs found in database")
            print()
            print(" This could mean:")
            print("   1. No jobs have been submitted yet")
            print("   2. Database was recently recreated")
            return

        print(f"\n Found {len(jobs)} recent jobs:")
        print()

        for i, job in enumerate(jobs, 1):
            job_id, sample, status, step, error, created, r1_path, r2_path = job

            print(f"Job #{i}: {sample}")
            print(f"  Status: {status}")
            print(f"  Current Step: {step or 'Not started'}")
            print(f"  Created: {created}")

            # Check FASTQ files
            if r1_path:
                r1_exists = os.path.exists(r1_path)
                r2_exists = os.path.exists(r2_path) if r2_path else False

                print(f"\n   Input Files:")
                print(f"     R1: {r1_path}")
                print(f"        {' EXISTS' if r1_exists else '❌ NOT FOUND'}")

                if r2_path:
                    print(f"     R2: {r2_path}")
                    print(f"        {'EXISTS' if r2_exists else '❌ NOT FOUND'}")

                if not r1_exists or not r2_exists:
                    print("\n    Input files are missing!")
                    print("     This could be why the job failed.")

            # Show error if exists
            if error:
                print(f"\n   Error Message:")
                error_lines = error.split('\n')
                for line in error_lines[:10]:  # Show first 10 lines
                    print(f"     {line}")
                if len(error_lines) > 10:
                    print(f"     ... ({len(error_lines) - 10} more lines)")

            # Check for output
            results_dir = Path("results") / job_id
            if results_dir.exists():
                print(f"\n  Results directory exists: {results_dir}")

                # Check for specific outputs
                output_dir = results_dir / "output"
                if output_dir.exists():
                    germline_vcf = output_dir / "Germline_VCF"
                    if germline_vcf.exists():
                        tsv_files = list(germline_vcf.glob("*.tsv"))
                        vcf_files = list(germline_vcf.glob("*.vcf*"))
                        print(f"     TSV files: {len(tsv_files)}")
                        print(f"     VCF files: {len(vcf_files)}")
            else:
                print(f"\n   No results directory found")

            print()
            print("-" * 60)
            print()

    except sqlite3.OperationalError as e:
        print(f" Database error: {e}")
        print()
        print(" This might mean:")
        print("   1. Database schema is outdated")
        print("   2. Database is corrupted")
        print()
        print("   Try restarting the backend to reinitialize the database.")

    finally:
        conn.close()

    # Check configuration
    print("\n Configuration Check:")
    print()

    try:
        from config import settings

        print(f"  Nextflow script: {settings.NEXTFLOW_SCRIPT}")
        nf_exists = os.path.exists(settings.NEXTFLOW_SCRIPT)
        print(f"    {'EXISTS' if nf_exists else '❌ NOT FOUND'}")

        print(f"\n  Reference genome: {settings.REFERENCE_GENOME}")
        ref_exists = os.path.exists(settings.REFERENCE_GENOME)
        print(f"    {'EXISTS' if ref_exists else '❌ NOT FOUND'}")

        print(f"\n  Upload directory: {settings.UPLOAD_DIR}")
        upload_dir = Path(settings.UPLOAD_DIR)
        if upload_dir.exists():
            file_count = len(list(upload_dir.glob("**/*")))
            print(f"     EXISTS ({file_count} items)")
        else:
            print(f"     NOT FOUND")

        print(f"\n  Results directory: {settings.RESULTS_DIR}")
        results_dir = Path(settings.RESULTS_DIR)
        if results_dir.exists():
            dir_count = len(list(results_dir.glob("*")))
            print(f"     EXISTS ({dir_count} job directories)")
        else:
            print(f"     NOT FOUND")

    except Exception as e:
        print(f" Error loading configuration: {e}")

    print()
    print("=" * 60)
    print()
    print(" Common Issues & Solutions:")
    print()
    print("1. 'Cannot find script file' error:")
    print("   → Check NEXTFLOW_SCRIPT path in .env")
    print("   → Should be: NEXTFLOW_SCRIPT=../main.nf")
    print()
    print("2. 'Quality Control' step failure:")
    print("   → Check if FASTQ files are valid")
    print("   → Check if fastp tool is installed")
    print("   → Look at work directory logs")
    print()
    print("3. Input files missing:")
    print("   → Files might not have uploaded correctly")
    print("   → Check uploads/ directory")
    print()

if __name__ == "__main__":
    diagnose_jobs()
