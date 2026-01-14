#!/usr/bin/env python3
"""
MinIO Security Configuration Script
Sets up IAM policies and access controls for WES pipeline
Created: 2026-01-14
"""

from minio import Minio
from minio.error import S3Error
import json
import sys

# MinIO connection settings
MINIO_ENDPOINT = "localhost:9000"
MINIO_ACCESS_KEY = "admin"
MINIO_SECRET_KEY = "WES2026SecureMinIO!"
MINIO_SECURE = False

# Bucket policies - restrictive by default
BUCKET_POLICIES = {
    "wes-raw-data": {
        "Version": "2012-10-17",
        "Statement": [
            {
                "Effect": "Allow",
                "Principal": {"AWS": ["*"]},
                "Action": ["s3:PutObject"],
                "Resource": ["arn:aws:s3:::wes-raw-data/*"],
                "Condition": {
                    "StringEquals": {
                        "s3:x-amz-server-side-encryption": "AES256"
                    }
                }
            },
            {
                "Effect": "Allow",
                "Principal": {"AWS": ["*"]},
                "Action": ["s3:GetObject"],
                "Resource": ["arn:aws:s3:::wes-raw-data/*"]
            }
        ]
    },
    "wes-intermediate": {
        "Version": "2012-10-17",
        "Statement": [
            {
                "Effect": "Allow",
                "Principal": {"AWS": ["*"]},
                "Action": ["s3:PutObject", "s3:GetObject", "s3:DeleteObject"],
                "Resource": ["arn:aws:s3:::wes-intermediate/*"]
            }
        ]
    },
    "wes-results": {
        "Version": "2012-10-17",
        "Statement": [
            {
                "Effect": "Allow",
                "Principal": {"AWS": ["*"]},
                "Action": ["s3:PutObject", "s3:GetObject"],
                "Resource": ["arn:aws:s3:::wes-results/*"]
            }
        ]
    },
    "wes-archives": {
        "Version": "2012-10-17",
        "Statement": [
            {
                "Effect": "Allow",
                "Principal": {"AWS": ["*"]},
                "Action": ["s3:PutObject", "s3:GetObject"],
                "Resource": ["arn:aws:s3:::wes-archives/*"]
            }
        ]
    },
    "wes-reference": {
        "Version": "2012-10-17",
        "Statement": [
            {
                "Effect": "Allow",
                "Principal": {"AWS": ["*"]},
                "Action": ["s3:GetObject"],
                "Resource": ["arn:aws:s3:::wes-reference/*"]
            }
        ]
    },
    "wes-logs": {
        "Version": "2012-10-17",
        "Statement": [
            {
                "Effect": "Allow",
                "Principal": {"AWS": ["*"]},
                "Action": ["s3:PutObject", "s3:GetObject"],
                "Resource": ["arn:aws:s3:::wes-logs/*"]
            }
        ]
    }
}

def configure_security():
    """
    Configure security policies for all MinIO buckets
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

        print("✓ Connected to MinIO\n")

        # Apply policies to each bucket
        configured_count = 0

        for bucket_name, policy in BUCKET_POLICIES.items():
            try:
                # Convert policy to JSON string
                policy_json = json.dumps(policy)

                # Set bucket policy
                client.set_bucket_policy(bucket_name, policy_json)

                print(f"✓ Configured security policy for '{bucket_name}'")

                # Display policy summary
                statements = policy.get("Statement", [])
                for stmt in statements:
                    actions = stmt.get("Action", [])
                    if isinstance(actions, str):
                        actions = [actions]
                    print(f"  - {stmt.get('Effect', 'N/A')}: {', '.join(actions)}")

                print()
                configured_count += 1

            except S3Error as e:
                print(f"✗ Error configuring '{bucket_name}': {e}")
                # Continue with other buckets
                continue

        # Summary
        print("=" * 60)
        print(f"Security Configuration Complete!")
        print(f"  Configured: {configured_count}/{len(BUCKET_POLICIES)}")
        print("=" * 60)

        # Security recommendations
        print("\n🔒 Security Recommendations:")
        print("  1. Change default admin credentials immediately")
        print("  2. Create service-specific users with limited permissions")
        print("  3. Enable TLS/SSL for production deployments")
        print("  4. Implement network firewall rules")
        print("  5. Enable audit logging")
        print("  6. Rotate credentials regularly")

        return 0 if configured_count == len(BUCKET_POLICIES) else 1

    except Exception as e:
        print(f"✗ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(configure_security())
