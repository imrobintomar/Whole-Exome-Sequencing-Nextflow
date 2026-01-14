# MinIO Object Storage Implementation Plan
## Whole Exome Sequencing SaaS Platform

**Date:** 2026-01-14
**Version:** 1.0
**Target Architecture:** MinIO-based Object Storage for WES Pipeline

---

## 📋 Table of Contents

1. [Executive Summary](#executive-summary)
2. [Current vs Target Architecture](#current-vs-target-architecture)
3. [MinIO Infrastructure Setup](#minio-infrastructure-setup)
4. [Storage Strategy & Buckets](#storage-strategy--buckets)
5. [Backend Integration](#backend-integration)
6. [Nextflow Pipeline Integration](#nextflow-pipeline-integration)
7. [Security & Access Control](#security--access-control)
8. [Data Lifecycle Management](#data-lifecycle-management)
9. [Performance Optimization](#performance-optimization)
10. [Monitoring & Logging](#monitoring--logging)
11. [Implementation Roadmap](#implementation-roadmap)
12. [Cost Estimation](#cost-estimation)

---

## 1. Executive Summary

### Why MinIO?

MinIO is the ideal choice for genomics SaaS platforms:

✅ **S3-Compatible API** - Works with AWS SDK, easy migration path
✅ **High Performance** - Optimized for large files (BAM/VCF/FASTQ)
✅ **Cost Effective** - No egress fees, predictable pricing
✅ **Self-Hosted** - Full data control, HIPAA-compliant deployment
✅ **Kubernetes Native** - Cloud-agnostic, scalable
✅ **Encryption** - At-rest and in-transit encryption
✅ **Versioning** - Object versioning for compliance

### Implementation Goals

1. **Replace local file storage** with MinIO object storage
2. **Implement bucket-based organization** (raw data, results, archives)
3. **Enable secure multi-tenant access** with IAM policies
4. **Automate data lifecycle** (hot → warm → cold → delete)
5. **Integrate with Nextflow** for seamless pipeline execution
6. **Ensure HIPAA compliance** with encryption and audit logging

---

## 2. Current vs Target Architecture

### Current Architecture (File-Based)

```
┌─────────────────────────────────────────────────────────────┐
│                    Local File System                         │
├─────────────────────────────────────────────────────────────┤
│  /media/drprabudh/m3/Nextflow-Script/WholeExome/backend/    │
│    ├── results/                                              │
│    │   └── {job_id}/                                        │
│    │       ├── input/  (FASTQ files)                        │
│    │       └── output/                                       │
│    │           ├── filtered_fastp/                           │
│    │           ├── Mapsam/ (BAM files)                      │
│    │           ├── Germline_VCF/ (VCF files)                │
│    │           └── Final_Annotated/ (TSV files)             │
│    └── work/  (Nextflow work directory)                     │
│                                                              │
│  Issues:                                                     │
│  ❌ Not scalable across multiple servers                    │
│  ❌ No redundancy or backup strategy                        │
│  ❌ Limited multi-tenancy support                           │
│  ❌ Difficult to implement tiered storage                   │
│  ❌ Manual cleanup required                                 │
└─────────────────────────────────────────────────────────────┘
```

### Target Architecture (MinIO-Based)

```
┌─────────────────────────────────────────────────────────────────────────┐
│                         MinIO Object Storage                             │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                          │
│  Buckets:                                                                │
│  ┌─────────────────────────────────────────────────────────────┐       │
│  │ wes-raw-data (HOT Storage - 7 days)                         │       │
│  │   └── {user_id}/{job_id}/                                   │       │
│  │       ├── {sample}_1.fastq.gz                               │       │
│  │       └── {sample}_2.fastq.gz                               │       │
│  └─────────────────────────────────────────────────────────────┘       │
│                                                                          │
│  ┌─────────────────────────────────────────────────────────────┐       │
│  │ wes-intermediate (WARM Storage - 30 days)                   │       │
│  │   └── {user_id}/{job_id}/                                   │       │
│  │       ├── filtered_fastp/                                    │       │
│  │       ├── bam/ (BAM files)                                  │       │
│  │       └── vcf/ (Raw VCF files)                              │       │
│  └─────────────────────────────────────────────────────────────┘       │
│                                                                          │
│  ┌─────────────────────────────────────────────────────────────┐       │
│  │ wes-results (HOT Storage - 90 days)                         │       │
│  │   └── {user_id}/{job_id}/                                   │       │
│  │       ├── {sample}_Final_.txt (annotated variants)          │       │
│  │       ├── {sample}_report.pdf (clinical report)             │       │
│  │       └── metadata.json                                      │       │
│  └─────────────────────────────────────────────────────────────┘       │
│                                                                          │
│  ┌─────────────────────────────────────────────────────────────┐       │
│  │ wes-archives (COLD Storage - 7 years HIPAA)                 │       │
│  │   └── {year}/{month}/{user_id}/{job_id}.tar.gz             │       │
│  └─────────────────────────────────────────────────────────────┘       │
│                                                                          │
│  ┌─────────────────────────────────────────────────────────────┐       │
│  │ wes-reference (Read-Only - No expiration)                   │       │
│  │   ├── hg38/                                                  │       │
│  │   │   ├── hg38.fa                                           │       │
│  │   │   ├── hg38.fa.fai                                       │       │
│  │   │   └── hg38.dict                                         │       │
│  │   ├── annovar_db/                                           │       │
│  │   └── known_sites/                                          │       │
│  └─────────────────────────────────────────────────────────────┘       │
│                                                                          │
│  ┌─────────────────────────────────────────────────────────────┐       │
│  │ wes-logs (WARM Storage - 1 year)                            │       │
│  │   └── {date}/{job_id}/                                      │       │
│  │       ├── trace.txt                                          │       │
│  │       ├── report.html                                        │       │
│  │       └── audit_logs.json                                    │       │
│  └─────────────────────────────────────────────────────────────┘       │
│                                                                          │
│  Features:                                                               │
│  ✅ Distributed, redundant storage (erasure coding)                    │
│  ✅ Multi-server scalability                                           │
│  ✅ Automated lifecycle policies (hot → warm → cold → delete)          │
│  ✅ Per-user IAM policies and access control                           │
│  ✅ Encryption at rest (AES-256) and in transit (TLS)                  │
│  ✅ Object versioning for compliance                                   │
│  ✅ S3-compatible API (boto3, AWS CLI, mc)                             │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## 3. MinIO Infrastructure Setup

### 3.1 Installation Options

#### Option A: Standalone MinIO (Development/Small Scale)

```bash
# Download MinIO server
wget https://dl.min.io/server/minio/release/linux-amd64/minio
chmod +x minio
sudo mv minio /usr/local/bin/

# Download MinIO client (mc)
wget https://dl.min.io/client/mc/release/linux-amd64/mc
chmod +x mc
sudo mv mc /usr/local/bin/

# Create data directory
sudo mkdir -p /mnt/minio/data
sudo chown -R $USER:$USER /mnt/minio

# Create systemd service
sudo tee /etc/systemd/system/minio.service << 'EOF'
[Unit]
Description=MinIO Object Storage
Documentation=https://docs.min.io
After=network-online.target
Wants=network-online.target

[Service]
Type=notify
WorkingDirectory=/usr/local/
User=minio
Group=minio

# Environment variables
Environment="MINIO_ROOT_USER=admin"
Environment="MINIO_ROOT_PASSWORD=YourStrongPassword123!"
Environment="MINIO_VOLUMES=/mnt/minio/data"
Environment="MINIO_OPTS=--console-address :9001"

# Enable encryption
Environment="MINIO_KMS_SECRET_KEY=my-minio-key:CHANGE-THIS-TO-A-SECURE-32-CHAR-KEY"

ExecStartPre=/bin/bash -c "if [ -z \"${MINIO_VOLUMES}\" ]; then echo \"Variable MINIO_VOLUMES not set\"; exit 1; fi"
ExecStart=/usr/local/bin/minio server $MINIO_OPTS $MINIO_VOLUMES

Restart=always
LimitNOFILE=65536
TasksMax=infinity
TimeoutStopSec=infinity
SendSIGKILL=no

[Install]
WantedBy=multi-user.target
EOF

# Create minio user
sudo useradd -r -s /sbin/nologin minio
sudo chown -R minio:minio /mnt/minio

# Start MinIO
sudo systemctl daemon-reload
sudo systemctl enable minio
sudo systemctl start minio

# Check status
sudo systemctl status minio
```

**Access Points:**
- API Endpoint: `http://localhost:9000`
- Console UI: `http://localhost:9001`

#### Option B: Distributed MinIO (Production/High Availability)

```bash
# 4-node distributed setup with erasure coding
# Deploy on 4 servers: node1, node2, node3, node4

# On each node, create systemd service with distributed volumes
Environment="MINIO_VOLUMES=http://node{1...4}/mnt/disk{1...4}/minio"

# Example for 4 nodes × 4 disks = 16 drives total
# Provides N/2 redundancy (can lose 8 drives)
```

**Production Configuration:**
```yaml
# docker-compose.yml for distributed MinIO with load balancer
version: '3.8'

services:
  minio1:
    image: minio/minio:latest
    hostname: minio1
    volumes:
      - /mnt/disk1:/data1
      - /mnt/disk2:/data2
    environment:
      MINIO_ROOT_USER: admin
      MINIO_ROOT_PASSWORD: ${MINIO_ROOT_PASSWORD}
      MINIO_VOLUMES: "http://minio{1...4}/data{1...2}"
      MINIO_OPTS: "--console-address :9001"
    command: server --address ":9000" --console-address ":9001"
    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:9000/minio/health/live"]
      interval: 30s
      timeout: 20s
      retries: 3

  minio2:
    # ... (similar configuration for minio2, minio3, minio4)

  nginx:
    image: nginx:alpine
    ports:
      - "9000:9000"
      - "9001:9001"
    volumes:
      - ./nginx.conf:/etc/nginx/nginx.conf:ro
    depends_on:
      - minio1
      - minio2
      - minio3
      - minio4
```

### 3.2 Initial Configuration

```bash
# Configure mc client
mc alias set myminio http://localhost:9000 admin YourStrongPassword123!

# Verify connection
mc admin info myminio

# Set bucket versioning (HIPAA compliance)
mc version enable myminio

# Configure encryption
mc admin kms key create myminio wes-encryption-key
```

---

## 4. Storage Strategy & Buckets

### 4.1 Bucket Structure

```bash
#!/bin/bash
# create_buckets.sh - Initialize MinIO buckets

MINIO_ALIAS="myminio"

# Create buckets
mc mb ${MINIO_ALIAS}/wes-raw-data
mc mb ${MINIO_ALIAS}/wes-intermediate
mc mb ${MINIO_ALIAS}/wes-results
mc mb ${MINIO_ALIAS}/wes-archives
mc mb ${MINIO_ALIAS}/wes-reference
mc mb ${MINIO_ALIAS}/wes-logs

# Enable versioning (compliance requirement)
mc version enable ${MINIO_ALIAS}/wes-raw-data
mc version enable ${MINIO_ALIAS}/wes-results
mc version enable ${MINIO_ALIAS}/wes-archives

# Set encryption
mc encrypt set sse-s3 ${MINIO_ALIAS}/wes-raw-data
mc encrypt set sse-s3 ${MINIO_ALIAS}/wes-intermediate
mc encrypt set sse-s3 ${MINIO_ALIAS}/wes-results
mc encrypt set sse-s3 ${MINIO_ALIAS}/wes-archives

# Configure lifecycle policies
# Raw data: Delete after 7 days
cat > raw-data-lifecycle.json << 'EOF'
{
  "Rules": [
    {
      "ID": "DeleteRawDataAfter7Days",
      "Status": "Enabled",
      "Expiration": {
        "Days": 7
      }
    }
  ]
}
EOF
mc ilm import ${MINIO_ALIAS}/wes-raw-data < raw-data-lifecycle.json

# Intermediate files: Delete after 30 days
cat > intermediate-lifecycle.json << 'EOF'
{
  "Rules": [
    {
      "ID": "DeleteIntermediateAfter30Days",
      "Status": "Enabled",
      "Expiration": {
        "Days": 30
      }
    }
  ]
}
EOF
mc ilm import ${MINIO_ALIAS}/wes-intermediate < intermediate-lifecycle.json

# Results: Transition to archives after 90 days, delete after 7 years
cat > results-lifecycle.json << 'EOF'
{
  "Rules": [
    {
      "ID": "ArchiveResultsAfter90Days",
      "Status": "Enabled",
      "Transition": {
        "Days": 90,
        "StorageClass": "GLACIER"
      }
    },
    {
      "ID": "DeleteResultsAfter7Years",
      "Status": "Enabled",
      "Expiration": {
        "Days": 2555
      }
    }
  ]
}
EOF
mc ilm import ${MINIO_ALIAS}/wes-results < results-lifecycle.json

echo "✓ Buckets created and configured"
```

### 4.2 Bucket Policies

```json
// wes-raw-data-policy.json
// Users can only read/write their own data
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Principal": {
        "AWS": ["*"]
      },
      "Action": [
        "s3:PutObject",
        "s3:GetObject",
        "s3:DeleteObject"
      ],
      "Resource": [
        "arn:aws:s3:::wes-raw-data/${aws:userid}/*"
      ],
      "Condition": {
        "StringLike": {
          "s3:prefix": ["${aws:userid}/*"]
        }
      }
    }
  ]
}
```

```bash
# Apply bucket policy
mc policy set-json wes-raw-data-policy.json myminio/wes-raw-data
```

---

## 5. Backend Integration

### 5.1 Python MinIO Client Setup

**Install Dependencies:**
```bash
pip install minio boto3 python-dotenv
```

**Environment Configuration:**
```bash
# .env
MINIO_ENDPOINT=localhost:9000
MINIO_ACCESS_KEY=admin
MINIO_SECRET_KEY=YourStrongPassword123!
MINIO_SECURE=False  # Set to True for HTTPS
MINIO_REGION=us-east-1

# Bucket names
MINIO_BUCKET_RAW=wes-raw-data
MINIO_BUCKET_INTERMEDIATE=wes-intermediate
MINIO_BUCKET_RESULTS=wes-results
MINIO_BUCKET_ARCHIVES=wes-archives
MINIO_BUCKET_REFERENCE=wes-reference
MINIO_BUCKET_LOGS=wes-logs
```

### 5.2 MinIO Service Layer

**Create `/backend/services/minio_service.py`:**

```python
"""
MinIO Object Storage Service
Handles all S3-compatible object storage operations
"""
from minio import Minio
from minio.error import S3Error
from minio.commonconfig import CopySource
from datetime import timedelta
from pathlib import Path
import os
import logging
from typing import Optional, List, BinaryIO
from config import settings

logger = logging.getLogger(__name__)


class MinIOService:
    """MinIO object storage service"""

    def __init__(self):
        """Initialize MinIO client"""
        self.client = Minio(
            settings.MINIO_ENDPOINT,
            access_key=settings.MINIO_ACCESS_KEY,
            secret_key=settings.MINIO_SECRET_KEY,
            secure=settings.MINIO_SECURE
        )

        # Verify connection
        try:
            self.client.list_buckets()
            logger.info("✓ MinIO connection established")
        except S3Error as e:
            logger.error(f"MinIO connection failed: {e}")
            raise

    def ensure_bucket_exists(self, bucket_name: str):
        """Create bucket if it doesn't exist"""
        try:
            if not self.client.bucket_exists(bucket_name):
                self.client.make_bucket(bucket_name)
                logger.info(f"Created bucket: {bucket_name}")
        except S3Error as e:
            logger.error(f"Failed to create bucket {bucket_name}: {e}")
            raise

    def upload_file(
        self,
        bucket_name: str,
        object_name: str,
        file_path: str,
        content_type: str = "application/octet-stream",
        metadata: dict = None
    ) -> bool:
        """
        Upload file to MinIO

        Args:
            bucket_name: Target bucket
            object_name: Object key/path in bucket
            file_path: Local file path
            content_type: MIME type
            metadata: Custom metadata dict

        Returns:
            bool: Upload success
        """
        try:
            self.ensure_bucket_exists(bucket_name)

            # Upload file
            self.client.fput_object(
                bucket_name,
                object_name,
                file_path,
                content_type=content_type,
                metadata=metadata
            )

            file_size = Path(file_path).stat().st_size
            logger.info(
                f"✓ Uploaded {file_path} → s3://{bucket_name}/{object_name} "
                f"({file_size / (1024**2):.2f} MB)"
            )
            return True

        except S3Error as e:
            logger.error(f"Upload failed: {e}")
            return False

    def upload_stream(
        self,
        bucket_name: str,
        object_name: str,
        data: BinaryIO,
        length: int,
        content_type: str = "application/octet-stream",
        metadata: dict = None
    ) -> bool:
        """Upload data from stream/file-like object"""
        try:
            self.ensure_bucket_exists(bucket_name)

            self.client.put_object(
                bucket_name,
                object_name,
                data,
                length,
                content_type=content_type,
                metadata=metadata
            )

            logger.info(f"✓ Uploaded stream → s3://{bucket_name}/{object_name}")
            return True

        except S3Error as e:
            logger.error(f"Stream upload failed: {e}")
            return False

    def download_file(
        self,
        bucket_name: str,
        object_name: str,
        file_path: str
    ) -> bool:
        """Download file from MinIO"""
        try:
            # Create parent directory
            Path(file_path).parent.mkdir(parents=True, exist_ok=True)

            # Download file
            self.client.fget_object(bucket_name, object_name, file_path)

            file_size = Path(file_path).stat().st_size
            logger.info(
                f"✓ Downloaded s3://{bucket_name}/{object_name} → {file_path} "
                f"({file_size / (1024**2):.2f} MB)"
            )
            return True

        except S3Error as e:
            logger.error(f"Download failed: {e}")
            return False

    def get_object(self, bucket_name: str, object_name: str):
        """Get object as stream"""
        try:
            return self.client.get_object(bucket_name, object_name)
        except S3Error as e:
            logger.error(f"Get object failed: {e}")
            return None

    def delete_object(self, bucket_name: str, object_name: str) -> bool:
        """Delete object from MinIO"""
        try:
            self.client.remove_object(bucket_name, object_name)
            logger.info(f"✓ Deleted s3://{bucket_name}/{object_name}")
            return True
        except S3Error as e:
            logger.error(f"Delete failed: {e}")
            return False

    def delete_objects(self, bucket_name: str, object_names: List[str]) -> bool:
        """Delete multiple objects (batch operation)"""
        try:
            from minio.deleteobjects import DeleteObject

            delete_objects = [DeleteObject(name) for name in object_names]
            errors = self.client.remove_objects(bucket_name, delete_objects)

            # Check for errors
            error_list = list(errors)
            if error_list:
                logger.error(f"Some deletions failed: {error_list}")
                return False

            logger.info(f"✓ Deleted {len(object_names)} objects from {bucket_name}")
            return True

        except S3Error as e:
            logger.error(f"Batch delete failed: {e}")
            return False

    def list_objects(
        self,
        bucket_name: str,
        prefix: str = "",
        recursive: bool = True
    ) -> List[str]:
        """List objects in bucket"""
        try:
            objects = self.client.list_objects(
                bucket_name,
                prefix=prefix,
                recursive=recursive
            )
            return [obj.object_name for obj in objects]
        except S3Error as e:
            logger.error(f"List objects failed: {e}")
            return []

    def get_presigned_url(
        self,
        bucket_name: str,
        object_name: str,
        expires: timedelta = timedelta(hours=1)
    ) -> Optional[str]:
        """
        Generate presigned URL for temporary download access

        Args:
            bucket_name: Bucket name
            object_name: Object key
            expires: URL expiration time (default 1 hour)

        Returns:
            str: Presigned URL or None
        """
        try:
            url = self.client.presigned_get_object(
                bucket_name,
                object_name,
                expires=expires
            )
            logger.info(f"Generated presigned URL for {object_name} (expires in {expires})")
            return url
        except S3Error as e:
            logger.error(f"Presigned URL generation failed: {e}")
            return None

    def get_upload_url(
        self,
        bucket_name: str,
        object_name: str,
        expires: timedelta = timedelta(hours=1)
    ) -> Optional[str]:
        """Generate presigned URL for direct upload"""
        try:
            url = self.client.presigned_put_object(
                bucket_name,
                object_name,
                expires=expires
            )
            return url
        except S3Error as e:
            logger.error(f"Upload URL generation failed: {e}")
            return None

    def copy_object(
        self,
        src_bucket: str,
        src_object: str,
        dst_bucket: str,
        dst_object: str
    ) -> bool:
        """Copy object between buckets"""
        try:
            self.client.copy_object(
                dst_bucket,
                dst_object,
                CopySource(src_bucket, src_object)
            )
            logger.info(f"✓ Copied {src_bucket}/{src_object} → {dst_bucket}/{dst_object}")
            return True
        except S3Error as e:
            logger.error(f"Copy failed: {e}")
            return False

    def get_object_stat(self, bucket_name: str, object_name: str) -> dict:
        """Get object metadata"""
        try:
            stat = self.client.stat_object(bucket_name, object_name)
            return {
                "size": stat.size,
                "last_modified": stat.last_modified,
                "etag": stat.etag,
                "content_type": stat.content_type,
                "metadata": stat.metadata
            }
        except S3Error as e:
            logger.error(f"Stat failed: {e}")
            return None

    def upload_directory(
        self,
        bucket_name: str,
        local_dir: str,
        prefix: str = ""
    ) -> bool:
        """Upload entire directory to MinIO"""
        try:
            local_path = Path(local_dir)
            if not local_path.exists():
                logger.error(f"Directory not found: {local_dir}")
                return False

            uploaded = 0
            for file_path in local_path.rglob("*"):
                if file_path.is_file():
                    # Calculate object key
                    relative_path = file_path.relative_to(local_path)
                    object_name = f"{prefix}/{relative_path}" if prefix else str(relative_path)

                    # Upload file
                    if self.upload_file(bucket_name, object_name, str(file_path)):
                        uploaded += 1

            logger.info(f"✓ Uploaded {uploaded} files from {local_dir} to {bucket_name}/{prefix}")
            return True

        except Exception as e:
            logger.error(f"Directory upload failed: {e}")
            return False

    def download_directory(
        self,
        bucket_name: str,
        prefix: str,
        local_dir: str
    ) -> bool:
        """Download directory from MinIO"""
        try:
            objects = self.list_objects(bucket_name, prefix=prefix)

            downloaded = 0
            for object_name in objects:
                # Calculate local path
                local_path = Path(local_dir) / object_name.replace(prefix, "", 1).lstrip("/")

                # Download file
                if self.download_file(bucket_name, object_name, str(local_path)):
                    downloaded += 1

            logger.info(f"✓ Downloaded {downloaded} files from {bucket_name}/{prefix} to {local_dir}")
            return True

        except Exception as e:
            logger.error(f"Directory download failed: {e}")
            return False


# Singleton instance
minio_service = MinIOService()
```

### 5.3 Update Backend Config

**Update `/backend/config.py`:**

```python
from pydantic_settings import BaseSettings
from typing import Optional

class Settings(BaseSettings):
    # ... existing settings ...

    # MinIO Configuration
    MINIO_ENDPOINT: str = "localhost:9000"
    MINIO_ACCESS_KEY: str = "admin"
    MINIO_SECRET_KEY: str = "YourStrongPassword123!"
    MINIO_SECURE: bool = False  # True for HTTPS
    MINIO_REGION: str = "us-east-1"

    # Bucket names
    MINIO_BUCKET_RAW: str = "wes-raw-data"
    MINIO_BUCKET_INTERMEDIATE: str = "wes-intermediate"
    MINIO_BUCKET_RESULTS: str = "wes-results"
    MINIO_BUCKET_ARCHIVES: str = "wes-archives"
    MINIO_BUCKET_REFERENCE: str = "wes-reference"
    MINIO_BUCKET_LOGS: str = "wes-logs"

    class Config:
        env_file = ".env"

settings = Settings()
```

### 5.4 Update Job Submission Endpoint

**Update `/backend/main.py`:**

```python
from fastapi import FastAPI, UploadFile, File, Depends
from services.minio_service import minio_service
from services.audit_service import AuditService
from datetime import datetime
import uuid

@app.post("/jobs/submit")
async def submit_wes_job(
    sample_name: str,
    r1_file: UploadFile = File(...),
    r2_file: UploadFile = File(...),
    user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """Submit WES analysis job with MinIO storage"""

    # Generate job ID
    job_id = str(uuid.uuid4())

    # Define object paths in MinIO
    user_prefix = f"{user.firebase_uid}/{job_id}"
    r1_object = f"{user_prefix}/{sample_name}_1.fastq.gz"
    r2_object = f"{user_prefix}/{sample_name}_2.fastq.gz"

    try:
        # Upload FASTQ files to MinIO
        logger.info(f"Uploading FASTQ files for job {job_id}")

        # Upload R1
        r1_success = minio_service.upload_stream(
            bucket_name=settings.MINIO_BUCKET_RAW,
            object_name=r1_object,
            data=r1_file.file,
            length=r1_file.size,
            content_type="application/gzip",
            metadata={
                "user_id": user.firebase_uid,
                "job_id": job_id,
                "sample_name": sample_name,
                "uploaded_at": datetime.utcnow().isoformat()
            }
        )

        # Upload R2
        r2_success = minio_service.upload_stream(
            bucket_name=settings.MINIO_BUCKET_RAW,
            object_name=r2_object,
            data=r2_file.file,
            length=r2_file.size,
            content_type="application/gzip",
            metadata={
                "user_id": user.firebase_uid,
                "job_id": job_id,
                "sample_name": sample_name,
                "uploaded_at": datetime.utcnow().isoformat()
            }
        )

        if not (r1_success and r2_success):
            raise Exception("Failed to upload FASTQ files to MinIO")

        # Create job in database
        job = Job(
            job_id=job_id,
            user_id=user.id,
            sample_name=sample_name,
            status=JobStatus.PENDING,
            fastq_r1_path=f"s3://{settings.MINIO_BUCKET_RAW}/{r1_object}",
            fastq_r2_path=f"s3://{settings.MINIO_BUCKET_RAW}/{r2_object}",
            created_at=datetime.utcnow()
        )
        db.add(job)
        db.commit()

        # Audit log
        AuditService.log(
            db=db,
            user_id=user.firebase_uid,
            user_email=user.email,
            action="JOB_SUBMITTED",
            resource_type="job",
            resource_id=job_id,
            details={
                "sample_name": sample_name,
                "r1_size_mb": r1_file.size / (1024**2),
                "r2_size_mb": r2_file.size / (1024**2),
                "storage": "minio"
            },
            phi_accessed=False
        )

        # Trigger pipeline execution (async)
        # background_tasks.add_task(run_pipeline, job_id, db)

        return {
            "job_id": job_id,
            "status": "submitted",
            "message": "Job submitted successfully. Pipeline will start shortly.",
            "fastq_files": {
                "r1": f"s3://{settings.MINIO_BUCKET_RAW}/{r1_object}",
                "r2": f"s3://{settings.MINIO_BUCKET_RAW}/{r2_object}"
            }
        }

    except Exception as e:
        logger.error(f"Job submission failed: {e}")
        raise HTTPException(500, f"Job submission failed: {str(e)}")


@app.get("/jobs/{job_id}/download/{file_type}")
async def download_result_file(
    job_id: str,
    file_type: str,  # "vcf", "tsv", "report", "bam"
    user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """Generate presigned URL for file download"""

    # Verify job ownership
    job = db.query(Job).filter(Job.job_id == job_id, Job.user_id == user.id).first()
    if not job:
        raise HTTPException(404, "Job not found")

    # Map file type to object path
    file_mapping = {
        "vcf": f"{user.firebase_uid}/{job_id}/vcf/{job.sample_name}.vcf.gz",
        "tsv": f"{user.firebase_uid}/{job_id}/{job.sample_name}_Final_.txt",
        "report": f"{user.firebase_uid}/{job_id}/{job.sample_name}_report.pdf",
        "bam": f"{user.firebase_uid}/{job_id}/bam/{job.sample_name}_recall.bam"
    }

    object_name = file_mapping.get(file_type)
    if not object_name:
        raise HTTPException(400, f"Invalid file type: {file_type}")

    # Determine bucket
    bucket = settings.MINIO_BUCKET_RESULTS if file_type in ["tsv", "report"] else settings.MINIO_BUCKET_INTERMEDIATE

    # Generate presigned URL (1 hour expiration)
    url = minio_service.get_presigned_url(
        bucket_name=bucket,
        object_name=object_name,
        expires=timedelta(hours=1)
    )

    if not url:
        raise HTTPException(500, "Failed to generate download URL")

    # Audit log
    AuditService.log(
        db=db,
        user_id=user.firebase_uid,
        user_email=user.email,
        action="FILE_DOWNLOADED",
        resource_type="file",
        resource_id=f"{job_id}/{file_type}",
        details={"file_type": file_type, "object_name": object_name},
        phi_accessed=True,
        data_exported=True
    )

    return {
        "download_url": url,
        "expires_in_seconds": 3600,
        "file_type": file_type
    }
```

---

## 6. Nextflow Pipeline Integration

### 6.1 Nextflow AWS CLI Plugin Configuration

**Update `nextflow.config`:**

```groovy
// Enable AWS S3 support
plugins {
    id 'nf-amazon'
}

// MinIO as S3-compatible endpoint
aws {
    client {
        endpoint = 'http://localhost:9000'
        s3PathStyleAccess = true  // Required for MinIO
        signerOverride = 'S3SignerType'  // Use S3v4 signing
    }
    accessKey = 'admin'
    secretKey = 'YourStrongPassword123!'
    region = 'us-east-1'
}

// Update params to use S3 paths
params {
    // Input from MinIO
    input_fastq_r1 = "s3://wes-raw-data/${params.user_id}/${params.job_id}/*_1.fastq.gz"
    input_fastq_r2 = "s3://wes-raw-data/${params.user_id}/${params.job_id}/*_2.fastq.gz"

    // Output to MinIO
    output_dir = "s3://wes-results/${params.user_id}/${params.job_id}"

    // Reference files from MinIO
    reference = "s3://wes-reference/hg38/hg38.fa"
    reference_index = "s3://wes-reference/hg38/hg38.fa.fai"
    reference_dict = "s3://wes-reference/hg38/hg38.dict"

    // ANNOVAR database
    annovar_db = "s3://wes-reference/annovar_db/hg38_humandb"
}

// Work directory in MinIO (optional, can use local for performance)
workDir = "s3://wes-intermediate/${params.user_id}/${params.job_id}/work"

// Publish results to MinIO
process {
    publishDir = [
        path: "${params.output_dir}",
        mode: 'copy',
        overwrite: false
    ]
}
```

### 6.2 Alternative: Hybrid Approach (Recommended)

For better performance, use **local work directory** and only store **final results** in MinIO:

```groovy
// Hybrid configuration
workDir = "/tmp/nextflow-work/${params.job_id}"  // Local for speed

process {
    // Stage input from MinIO to local
    stageInMode = 'copy'

    // Different publish strategies per process
    withName: 'fastpQC' {
        publishDir = [
            path: "s3://wes-intermediate/${params.user_id}/${params.job_id}/qc",
            mode: 'copy'
        ]
    }

    withName: 'haplotypeCaller' {
        publishDir = [
            path: "s3://wes-intermediate/${params.user_id}/${params.job_id}/vcf",
            mode: 'copy'
        ]
    }

    withName: 'addUniqueID' {
        publishDir = [
            path: "s3://wes-results/${params.user_id}/${params.job_id}",
            mode: 'copy'
        ]
    }
}
```

### 6.3 Update Main Workflow

**Update `main.nf`:**

```groovy
#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ... (existing includes and validation) ...

workflow {
    log.info """
    ╔════════════════════════════════════════╗
    ║   Whole Exome Sequencing Pipeline      ║
    ║        MinIO Object Storage            ║
    ╚════════════════════════════════════════╝

    Job ID:              ${params.job_id}
    User ID:             ${params.user_id}
    Sample:              ${params.sample_name}

    Storage Backend:     MinIO (S3-compatible)
    Input Bucket:        wes-raw-data
    Output Bucket:       wes-results
    Reference Bucket:    wes-reference

    Pipeline started at: ${new Date()}
    """

    // Download FASTQ files from MinIO to local work directory
    // (Nextflow handles this automatically with S3 paths)

    // Read paired-end FASTQ files from MinIO
    Channel
        .fromPath(params.input_fastq_r1)
        .map { r1 ->
            // Find corresponding R2 file
            def r2 = file(r1.toString().replace("_1.fastq.gz", "_2.fastq.gz"))
            def sample_id = r1.simpleName.replaceAll(/_1$/, '')
            tuple(sample_id, r1, r2)
        }
        .set { read_pairs }

    // Run pipeline (same as before)
    fastp_out = fastpQC(read_pairs)
    bwa_out = bwaMem(fastp_out)
    bam_sorted = sortSam(bwa_out)
    flagstat(bam_sorted)
    marked_bam = markDuplicates(bam_sorted)
    marked_sorted = sortSamPostDedup(marked_bam)
    recal_table = baseRecalibrator(marked_sorted)
    final_bam = applyBQSR(marked_sorted, recal_table)
    vcf_raw = haplotypeCaller(final_bam)
    annovar_txt = annovarAnnotate(vcf_raw)
    final_output = addUniqueID(annovar_txt)

    // Upload completion notification
    final_output.subscribe { sample_id, output_file ->
        log.info """
        ✓ Pipeline completed successfully at ${new Date()}

        Results uploaded to:
        s3://wes-results/${params.user_id}/${params.job_id}/${sample_id}_Final_.txt
        """
    }
}
```

### 6.4 Pipeline Execution Script

**Update `/backend/services/pipeline.py`:**

```python
import subprocess
import logging
from pathlib import Path
from services.minio_service import minio_service
from config import settings

logger = logging.getLogger(__name__)


class PipelineRunner:

    def __init__(self, job_id: str, user_id: str, sample_name: str):
        self.job_id = job_id
        self.user_id = user_id
        self.sample_name = sample_name

    def run_pipeline(self):
        """Execute Nextflow pipeline with MinIO backend"""

        # Nextflow command with MinIO parameters
        cmd = [
            "nextflow", "run", "main.nf",
            "--job_id", self.job_id,
            "--user_id", self.user_id,
            "--sample_name", self.sample_name,
            "-profile", "standard",
            "-with-trace", f"s3://{settings.MINIO_BUCKET_LOGS}/{self.job_id}/trace.txt",
            "-with-report", f"s3://{settings.MINIO_BUCKET_LOGS}/{self.job_id}/report.html",
            "-with-timeline", f"s3://{settings.MINIO_BUCKET_LOGS}/{self.job_id}/timeline.html",
            "-resume"
        ]

        # Set AWS credentials for MinIO
        env = os.environ.copy()
        env.update({
            "AWS_ACCESS_KEY_ID": settings.MINIO_ACCESS_KEY,
            "AWS_SECRET_ACCESS_KEY": settings.MINIO_SECRET_KEY,
            "AWS_REGION": settings.MINIO_REGION,
        })

        try:
            logger.info(f"Starting pipeline for job {self.job_id}")

            # Execute Nextflow
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                env=env,
                cwd=Path(__file__).parent.parent.parent  # Project root
            )

            # Stream output
            for line in process.stdout:
                logger.info(line.strip())

            # Wait for completion
            returncode = process.wait()

            if returncode == 0:
                logger.info(f"✓ Pipeline completed successfully for job {self.job_id}")
                return True
            else:
                logger.error(f"Pipeline failed with exit code {returncode}")
                return False

        except Exception as e:
            logger.error(f"Pipeline execution error: {e}")
            return False
```

---

## 7. Security & Access Control

### 7.1 User-Specific IAM Policies

**Create per-user MinIO policies:**

```bash
# create_user_policy.sh
USER_ID=$1

cat > user-${USER_ID}-policy.json << EOF
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Action": [
        "s3:ListBucket"
      ],
      "Resource": [
        "arn:aws:s3:::wes-raw-data",
        "arn:aws:s3:::wes-results"
      ],
      "Condition": {
        "StringLike": {
          "s3:prefix": ["${USER_ID}/*"]
        }
      }
    },
    {
      "Effect": "Allow",
      "Action": [
        "s3:GetObject",
        "s3:PutObject",
        "s3:DeleteObject"
      ],
      "Resource": [
        "arn:aws:s3:::wes-raw-data/${USER_ID}/*",
        "arn:aws:s3:::wes-results/${USER_ID}/*"
      ]
    }
  ]
}
EOF

# Create MinIO user
mc admin user add myminio user-${USER_ID} GeneratedPassword123!

# Attach policy
mc admin policy create myminio user-${USER_ID}-policy user-${USER_ID}-policy.json
mc admin policy attach myminio user-${USER_ID}-policy --user=user-${USER_ID}

echo "✓ Created user: user-${USER_ID}"
```

### 7.2 TLS/HTTPS Configuration

**Generate TLS certificates:**

```bash
# Generate self-signed certificate (for development)
openssl req -new -x509 -days 3650 -nodes \
    -out /etc/minio/certs/public.crt \
    -keyout /etc/minio/certs/private.key \
    -subj "/C=US/ST=State/L=City/O=Org/CN=minio.yourdomain.com"

# Update MinIO service to use HTTPS
# Update MINIO_OPTS in systemd service:
# MINIO_OPTS="--certs-dir /etc/minio/certs --console-address :9001"

# Restart MinIO
sudo systemctl restart minio
```

**Update backend config:**
```python
MINIO_ENDPOINT = "minio.yourdomain.com:9000"
MINIO_SECURE = True
```

### 7.3 Encryption at Rest

MinIO automatically encrypts objects when configured:

```bash
# Set default encryption for all buckets
mc encrypt set sse-s3 myminio/wes-raw-data
mc encrypt set sse-s3 myminio/wes-intermediate
mc encrypt set sse-s3 myminio/wes-results
mc encrypt set sse-s3 myminio/wes-archives

# Verify encryption
mc encrypt info myminio/wes-raw-data
```

---

## 8. Data Lifecycle Management

### 8.1 Automated Lifecycle Policies

**Already configured in section 4.1, but here's monitoring:**

```bash
# Monitor lifecycle transitions
mc admin trace myminio --verbose

# Check lifecycle rules
mc ilm ls myminio/wes-raw-data
mc ilm ls myminio/wes-intermediate
mc ilm ls myminio/wes-results
```

### 8.2 Manual Archival Script

```python
# archive_job.py
from services.minio_service import minio_service
from config import settings
import tarfile
from datetime import datetime
import tempfile

def archive_job(user_id: str, job_id: str):
    """Archive completed job to cold storage"""

    # Download all job files to temp directory
    with tempfile.TemporaryDirectory() as tmpdir:
        # Download intermediate files
        minio_service.download_directory(
            bucket_name=settings.MINIO_BUCKET_INTERMEDIATE,
            prefix=f"{user_id}/{job_id}",
            local_dir=f"{tmpdir}/intermediate"
        )

        # Download results
        minio_service.download_directory(
            bucket_name=settings.MINIO_BUCKET_RESULTS,
            prefix=f"{user_id}/{job_id}",
            local_dir=f"{tmpdir}/results"
        )

        # Create tar.gz archive
        archive_name = f"{job_id}.tar.gz"
        archive_path = f"{tmpdir}/{archive_name}"

        with tarfile.open(archive_path, "w:gz") as tar:
            tar.add(f"{tmpdir}/intermediate", arcname="intermediate")
            tar.add(f"{tmpdir}/results", arcname="results")

        # Upload archive
        year = datetime.now().year
        month = datetime.now().month

        minio_service.upload_file(
            bucket_name=settings.MINIO_BUCKET_ARCHIVES,
            object_name=f"{year}/{month:02d}/{user_id}/{archive_name}",
            file_path=archive_path,
            metadata={
                "user_id": user_id,
                "job_id": job_id,
                "archived_at": datetime.utcnow().isoformat()
            }
        )

        # Delete intermediate files (keep results for 90 days per lifecycle)
        objects = minio_service.list_objects(
            settings.MINIO_BUCKET_INTERMEDIATE,
            prefix=f"{user_id}/{job_id}"
        )
        minio_service.delete_objects(settings.MINIO_BUCKET_INTERMEDIATE, objects)

        logger.info(f"✓ Archived job {job_id} to cold storage")
```

---

## 9. Performance Optimization

### 9.1 Parallel Upload/Download

```python
from concurrent.futures import ThreadPoolExecutor, as_completed

def parallel_upload_files(file_list: List[tuple], bucket_name: str, max_workers: int = 10):
    """Upload multiple files in parallel"""

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = []

        for local_path, object_name in file_list:
            future = executor.submit(
                minio_service.upload_file,
                bucket_name,
                object_name,
                local_path
            )
            futures.append(future)

        # Wait for completion
        for future in as_completed(futures):
            try:
                result = future.result()
                if not result:
                    logger.error("Upload failed")
            except Exception as e:
                logger.error(f"Upload exception: {e}")
```

### 9.2 MinIO Server Tuning

```bash
# Increase file descriptor limits
ulimit -n 65536

# Update /etc/systemd/system/minio.service
[Service]
LimitNOFILE=65536
```

### 9.3 Caching Reference Files Locally

```groovy
// nextflow.config - Cache reference files locally
params {
    // Download reference from MinIO once, cache locally
    reference_cache_dir = "/opt/wes-cache/reference"

    reference = file("${params.reference_cache_dir}/hg38.fa").exists() ?
        "${params.reference_cache_dir}/hg38.fa" :
        "s3://wes-reference/hg38/hg38.fa"
}

// Pre-download script
process downloadReference {
    storeDir "${params.reference_cache_dir}"

    output:
    path "hg38.fa"

    script:
    """
    aws s3 cp s3://wes-reference/hg38/hg38.fa hg38.fa
    """
}
```

---

## 10. Monitoring & Logging

### 10.1 MinIO Prometheus Metrics

```bash
# Enable Prometheus metrics
mc admin prometheus generate myminio > prometheus.yml

# Add to prometheus.yml
scrape_configs:
  - job_name: 'minio'
    metrics_path: /minio/v2/metrics/cluster
    static_configs:
      - targets: ['localhost:9000']
```

### 10.2 Storage Usage Dashboard

```python
# Add to /backend/routes/monitoring.py

@router.get("/storage/usage")
async def get_storage_usage(
    admin_user: User = Depends(require_admin)
):
    """Get MinIO storage usage statistics"""

    from minio.helpers import get_target_url
    import requests

    # Get bucket sizes
    buckets = minio_service.client.list_buckets()

    usage = {}
    for bucket in buckets:
        objects = minio_service.list_objects(bucket.name, recursive=True)

        total_size = 0
        object_count = 0

        for obj_name in objects:
            stat = minio_service.get_object_stat(bucket.name, obj_name)
            if stat:
                total_size += stat['size']
                object_count += 1

        usage[bucket.name] = {
            "object_count": object_count,
            "total_size_gb": round(total_size / (1024**3), 2),
            "created_at": bucket.creation_date.isoformat()
        }

    return usage
```

---

## 11. Implementation Roadmap

### Phase 1: Infrastructure Setup (Week 1)

**Day 1-2: MinIO Installation**
- [ ] Install MinIO server (standalone or distributed)
- [ ] Configure systemd service
- [ ] Set up TLS/HTTPS
- [ ] Install MinIO client (mc)
- [ ] Configure backups

**Day 3-4: Bucket Configuration**
- [ ] Create all buckets (raw, intermediate, results, archives, reference, logs)
- [ ] Configure encryption (SSE-S3)
- [ ] Set up versioning
- [ ] Apply lifecycle policies
- [ ] Test bucket policies

**Day 5-7: Reference Data Migration**
- [ ] Upload hg38 reference genome
- [ ] Upload BWA indices
- [ ] Upload ANNOVAR databases
- [ ] Upload known sites VCFs
- [ ] Verify MD5 checksums

### Phase 2: Backend Integration (Week 2)

**Day 1-3: Python Integration**
- [ ] Install minio Python SDK
- [ ] Create MinIOService class
- [ ] Update config.py with MinIO settings
- [ ] Write unit tests for MinIO operations
- [ ] Test upload/download functionality

**Day 4-5: API Endpoints**
- [ ] Update job submission endpoint
- [ ] Implement presigned URL generation
- [ ] Add file download endpoints
- [ ] Update audit logging
- [ ] Test end-to-end workflow

**Day 6-7: Job Management**
- [ ] Update PipelineRunner for MinIO
- [ ] Implement job cleanup service
- [ ] Add archival functionality
- [ ] Test resume/retry logic

### Phase 3: Nextflow Integration (Week 3)

**Day 1-2: Nextflow Configuration**
- [ ] Install nf-amazon plugin
- [ ] Configure AWS settings for MinIO
- [ ] Update params for S3 paths
- [ ] Test S3 input/output

**Day 3-4: Process Updates**
- [ ] Update publishDir configurations
- [ ] Test each process with MinIO
- [ ] Optimize stageInMode/stageOutMode
- [ ] Handle large file transfers

**Day 5-7: Pipeline Testing**
- [ ] Run full pipeline with MinIO
- [ ] Test resume functionality
- [ ] Verify output correctness
- [ ] Performance benchmarking

### Phase 4: Security & Compliance (Week 4)

**Day 1-2: Access Control**
- [ ] Create IAM policies for users
- [ ] Set up bucket policies
- [ ] Test multi-tenant isolation
- [ ] Implement user-specific credentials

**Day 3-4: Encryption & Audit**
- [ ] Verify encryption at rest
- [ ] Configure TLS for all connections
- [ ] Set up audit logging
- [ ] Test compliance features

**Day 5-7: Monitoring & Alerts**
- [ ] Configure Prometheus metrics
- [ ] Set up Grafana dashboards
- [ ] Configure alerts
- [ ] Document operational procedures

### Phase 5: Production Deployment (Week 5-6)

**Week 5: Migration & Testing**
- [ ] Migrate existing data to MinIO
- [ ] Run parallel testing (old vs new)
- [ ] Performance tuning
- [ ] Load testing
- [ ] Disaster recovery testing

**Week 6: Go-Live**
- [ ] Final security audit
- [ ] Update documentation
- [ ] Train operations team
- [ ] Gradual rollout
- [ ] Monitor and optimize

---

## 12. Cost Estimation

### Hardware Requirements (Self-Hosted MinIO)

**Option A: Standalone (Dev/Small Scale)**
- **Storage:** 10 TB SSD ($1,500)
- **Server:** 32 GB RAM, 16 cores ($2,000)
- **Total:** ~$3,500 one-time

**Option B: Distributed 4-Node (Production)**
- **Storage:** 4 servers × 20 TB = 80 TB ($12,000)
- **Servers:** 4 × (64 GB RAM, 32 cores) ($16,000)
- **Network:** 10 Gbps switches ($2,000)
- **Total:** ~$30,000 one-time

### Operating Costs (Monthly)

| Item | Cost |
|------|------|
| Electricity (4 servers @ 500W) | $150 |
| Bandwidth (5 TB egress) | $0 (self-hosted) |
| Maintenance | $200 |
| **Total Monthly** | **$350** |

### Cost Comparison: MinIO vs AWS S3

**Assumptions:**
- 100 WES jobs/month
- 50 GB per job (raw FASTQ + intermediate + results)
- 5 TB total storage/month
- 10 TB egress/month

| Service | Storage Cost | Egress Cost | Total/Month |
|---------|--------------|-------------|-------------|
| **AWS S3 Standard** | $115 (5 TB) | $920 (10 TB) | **$1,035** |
| **AWS S3 Glacier** | $20 (5 TB) | $920 (10 TB) | **$940** |
| **MinIO Self-Hosted** | $0 (capex) | $0 | **$350** |

**Annual Savings with MinIO:** ~$8,220

**ROI:** Self-hosted MinIO pays for itself in 4-6 months

---

## 13. Backup & Disaster Recovery

### 13.1 MinIO Bucket Replication

```bash
# Set up replication to remote MinIO (DR site)
mc admin bucket remote add myminio/wes-results \
    https://dr-minio.example.com/wes-results-replica \
    --service replication \
    --access-key REMOTE_ACCESS_KEY \
    --secret-key REMOTE_SECRET_KEY

# Enable replication
mc replicate add myminio/wes-results \
    --remote-bucket wes-results-replica \
    --replicate "delete,delete-marker"
```

### 13.2 Automated Snapshots

```bash
# Daily backup script
#!/bin/bash
# /opt/scripts/minio-backup.sh

DATE=$(date +%Y%m%d)
BACKUP_DIR="/mnt/backups/minio/${DATE}"

# Create backup directory
mkdir -p ${BACKUP_DIR}

# Mirror all buckets
mc mirror myminio/wes-results ${BACKUP_DIR}/wes-results
mc mirror myminio/wes-archives ${BACKUP_DIR}/wes-archives

# Compress
tar -czf ${BACKUP_DIR}.tar.gz ${BACKUP_DIR}

# Remove old backups (keep 30 days)
find /mnt/backups/minio -type f -name "*.tar.gz" -mtime +30 -delete

echo "✓ Backup completed: ${BACKUP_DIR}.tar.gz"
```

**Add to cron:**
```bash
# crontab -e
0 2 * * * /opt/scripts/minio-backup.sh >> /var/log/minio-backup.log 2>&1
```

---

## 14. Testing Checklist

### Pre-Deployment Tests

- [ ] Upload 10 GB test file → verify speed and integrity
- [ ] Download 10 GB test file → verify speed and integrity
- [ ] Test parallel uploads (10 files simultaneously)
- [ ] Test presigned URL generation and expiration
- [ ] Verify encryption at rest (check etag format)
- [ ] Test TLS/HTTPS connections
- [ ] Verify bucket policies (try unauthorized access)
- [ ] Test lifecycle policies (manual date manipulation)
- [ ] Run full WES pipeline end-to-end
- [ ] Test pipeline resume after failure
- [ ] Verify audit logging
- [ ] Test disaster recovery (restore from backup)
- [ ] Load test (100 concurrent uploads)
- [ ] Monitor resource usage (CPU, RAM, disk I/O)

### Post-Deployment Monitoring

- [ ] Monitor daily storage growth
- [ ] Track upload/download latencies
- [ ] Check error rates
- [ ] Verify lifecycle transitions
- [ ] Review audit logs weekly
- [ ] Test disaster recovery quarterly

---

## 15. Documentation & Training

### Operational Documentation

**Create:**
- [ ] MinIO administration guide
- [ ] Backup and restore procedures
- [ ] Troubleshooting guide
- [ ] Security incident response plan
- [ ] Capacity planning guidelines

### Developer Documentation

**Update:**
- [ ] API documentation (presigned URLs, uploads)
- [ ] Nextflow configuration guide
- [ ] MinIO integration examples
- [ ] Testing procedures

---

## 16. Success Criteria

✅ **Infrastructure:**
- MinIO server running with 99.9% uptime
- All buckets created and configured
- TLS/HTTPS enabled
- Backup strategy implemented

✅ **Performance:**
- Upload speed: >100 MB/s for large files
- Download speed: >100 MB/s
- API response time: <200ms
- Pipeline execution time: No degradation vs local storage

✅ **Security:**
- Encryption at rest verified
- TLS in transit verified
- IAM policies tested
- Audit logging functional

✅ **Compliance:**
- HIPAA-compliant encryption
- Data retention policies active
- Audit logs retention: 1 year
- Disaster recovery tested

✅ **Integration:**
- Backend API fully functional
- Nextflow pipeline working end-to-end
- User upload/download working
- Monitoring dashboards deployed

---

## 17. Rollback Plan

If issues arise during deployment:

1. **Keep local storage active** during transition period
2. **Dual-write strategy:** Write to both local and MinIO initially
3. **Gradual migration:** Move 10% of users at a time
4. **Monitoring:** Track error rates, performance metrics
5. **Rollback trigger:** >5% error rate or >2x latency increase

**Rollback steps:**
```bash
# 1. Stop new jobs from using MinIO
export USE_MINIO=false

# 2. Re-route API to local storage
# Update config.py: USE_LOCAL_STORAGE=True

# 3. Keep MinIO data for forensics
# Don't delete buckets

# 4. Analyze logs and fix issues

# 5. Re-test and retry deployment
```

---

## 18. Conclusion

MinIO provides a **production-ready, cost-effective, HIPAA-compliant** object storage solution for your WES SaaS platform. This implementation plan ensures:

- **Scalability:** Handle growing data volumes
- **Performance:** High-speed uploads/downloads
- **Security:** Encryption, access control, audit logging
- **Compliance:** HIPAA, data retention, versioning
- **Cost Savings:** 8x cheaper than AWS S3
- **Flexibility:** S3-compatible API, easy migration

**Estimated Implementation Time:** 5-6 weeks
**Estimated Cost:** $3,500 - $30,000 (one-time hardware)
**Annual Savings:** ~$8,220 vs AWS S3

---

**Document Version:** 1.0
**Last Updated:** 2026-01-14
**Next Review:** After Phase 1 completion
**Owner:** DevOps & Platform Team
