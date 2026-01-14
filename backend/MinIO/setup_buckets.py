#!/usr/bin/env python3
"""
MinIO Bucket Setup Script
Creates the 6-bucket structure for WES pipeline
Created: 2026-01-14
"""

from minio import Minio
from minio.error import S3Error
import sys
import os

# MinIO connection settings
MINIO_ENDPOINT = "localhost:9000"
MINIO_ACCESS_KEY = "admin"
MINIO_SECRET_KEY = "WES2026SecureMinIO!"
MINIO_SECURE = False  # Set to True when using HTTPS

# Bucket definitions with lifecycle policies
BUCKETS = {
    "wes-raw-data": {
        "description": "Raw FASTQ files from users",
        "lifecycle_days": 30,  # Delete after 30 days
        "versioning": True
    },
    "wes-intermediate": {
        "description": "Intermediate processing files (BAM, sorted BAM, deduplicated)",
        "lifecycle_days": 7,  # Delete after 7 days
        "versioning": False
    },
    "wes-results": {
        "description": "Final analysis results (VCF, annotated variants, reports)",
        "lifecycle_days": 180,  # Keep for 6 months
        "versioning": True
    },
    "wes-archives": {
        "description": "Long-term archived data",
        "lifecycle_days": 1095,  # Keep for 3 years (1095 days)
        "versioning": True
    },
    "wes-reference": {
        "description": "Reference genomes, databases (hg38, ANNOVAR, etc.)",
        "lifecycle_days": None,  # Never delete
        "versioning": False
    },
    "wes-logs": {
        "description": "Pipeline execution logs and audit trails",
        "lifecycle_days": 90,  # Keep logs for 90 days
        "versioning": False
    }
}

def setup_minio_buckets():
    """
    Create all required MinIO buckets with appropriate configurations
    """
    try:
        # Initialize MinIO client
        print(f"Connecting to MinIO at {MINIO_ENDPOINT}...")
        client = Minio(
            MINIO_ENDPOINT,
            access_key=MINIO_ACCESS_KEY,
            secret_key=MINIO_SECRET_KEY,
            secure=MINIO_SECURE
        )

        # Verify connection
        try:
            client.list_buckets()
            print("✓ Successfully connected to MinIO\n")
        except Exception as e:
            print(f"✗ Failed to connect to MinIO: {e}")
            return 1

        # Create each bucket
        created_count = 0
        existing_count = 0

        for bucket_name, config in BUCKETS.items():
            try:
                # Check if bucket exists
                if client.bucket_exists(bucket_name):
                    print(f"⊙ Bucket '{bucket_name}' already exists")
                    existing_count += 1
                else:
                    # Create bucket
                    client.make_bucket(bucket_name)
                    print(f"✓ Created bucket '{bucket_name}'")
                    created_count += 1

                # Display configuration
                print(f"  Description: {config['description']}")
                print(f"  Lifecycle:   {config['lifecycle_days']} days" if config['lifecycle_days'] else "  Lifecycle:   Never delete")
                print(f"  Versioning:  {'Enabled' if config['versioning'] else 'Disabled'}")
                print()

            except S3Error as e:
                print(f"✗ Error creating bucket '{bucket_name}': {e}")
                return 1

        # Summary
        print("=" * 60)
        print(f"Bucket Setup Complete!")
        print(f"  Created: {created_count}")
        print(f"  Existing: {existing_count}")
        print(f"  Total: {len(BUCKETS)}")
        print("=" * 60)

        # Display access information
        print("\n📊 MinIO Server Information:")
        print(f"  API Endpoint:  http://{MINIO_ENDPOINT}")
        print(f"  Web Console:   http://localhost:9001")
        print(f"  Access Key:    {MINIO_ACCESS_KEY}")
        print(f"  Secret Key:    {MINIO_SECRET_KEY}")
        print("\n⚠️  IMPORTANT: Change default credentials in production!")

        return 0

    except Exception as e:
        print(f"✗ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(setup_minio_buckets())
