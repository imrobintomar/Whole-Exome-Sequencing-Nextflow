"""
MinIO Service Layer for WES Pipeline FastAPI Backend
Provides high-level operations for object storage management
Created: 2026-01-14
Updated: 2026-01-17 - Integrated with centralized config
"""

from minio import Minio
from minio.error import S3Error
from datetime import timedelta, datetime
from pathlib import Path
from typing import Optional, List, Dict, Tuple, BinaryIO
import logging
import os
import io

# Import centralized settings
try:
    from config import settings
    USE_CENTRALIZED_CONFIG = True
except ImportError:
    USE_CENTRALIZED_CONFIG = False

logger = logging.getLogger(__name__)


class MinIOConfig:
    """MinIO configuration settings - uses centralized config when available"""

    def __init__(self):
        if USE_CENTRALIZED_CONFIG:
            # Use centralized settings from config.py
            self.ENDPOINT = settings.MINIO_ENDPOINT
            self.ACCESS_KEY = settings.MINIO_ROOT_USER
            self.SECRET_KEY = settings.MINIO_ROOT_PASSWORD
            self.SECURE = settings.MINIO_SECURE
            self.ENABLED = settings.MINIO_ENABLED

            # Bucket names from config
            self.BUCKET_RAW_DATA = settings.MINIO_BUCKET_RAW_DATA
            self.BUCKET_INTERMEDIATE = settings.MINIO_BUCKET_INTERMEDIATE
            self.BUCKET_RESULTS = settings.MINIO_BUCKET_RESULTS
            self.BUCKET_ARCHIVES = settings.MINIO_BUCKET_ARCHIVES
            self.BUCKET_REFERENCE = settings.MINIO_BUCKET_REFERENCE
            self.BUCKET_LOGS = settings.MINIO_BUCKET_LOGS
        else:
            # Fallback to environment variables
            self.ENDPOINT = os.getenv("MINIO_ENDPOINT", "localhost:9000")
            self.ACCESS_KEY = os.getenv("MINIO_ROOT_USER", "admin")
            self.SECRET_KEY = os.getenv("MINIO_ROOT_PASSWORD", "WES2026SecureMinIO!")
            self.SECURE = os.getenv("MINIO_SECURE", "false").lower() == "true"
            self.ENABLED = os.getenv("MINIO_ENABLED", "false").lower() == "true"

            # Default bucket names
            self.BUCKET_RAW_DATA = "wes-raw-data"
            self.BUCKET_INTERMEDIATE = "wes-intermediate"
            self.BUCKET_RESULTS = "wes-results"
            self.BUCKET_ARCHIVES = "wes-archives"
            self.BUCKET_REFERENCE = "wes-reference"
            self.BUCKET_LOGS = "wes-logs"

    @property
    def all_buckets(self) -> List[str]:
        """Get list of all bucket names"""
        return [
            self.BUCKET_RAW_DATA,
            self.BUCKET_INTERMEDIATE,
            self.BUCKET_RESULTS,
            self.BUCKET_ARCHIVES,
            self.BUCKET_REFERENCE,
            self.BUCKET_LOGS
        ]


class MinIOService:
    """
    MinIO service layer for WES pipeline
    Handles all object storage operations
    """

    def __init__(self, config: MinIOConfig = None):
        """
        Initialize MinIO client

        Args:
            config: MinIOConfig instance (uses default if None)
        """
        self.config = config or MinIOConfig()

        try:
            self.client = Minio(
                self.config.ENDPOINT,
                access_key=self.config.ACCESS_KEY,
                secret_key=self.config.SECRET_KEY,
                secure=self.config.SECURE
            )
            logger.info(f"Connected to MinIO at {self.config.ENDPOINT}")
        except Exception as e:
            logger.error(f"Failed to connect to MinIO: {e}")
            raise

    # ==================== Upload Operations ====================

    def upload_file(
        self,
        bucket_name: str,
        object_name: str,
        file_path: str,
        metadata: Optional[Dict[str, str]] = None,
        content_type: str = "application/octet-stream"
    ) -> bool:
        """
        Upload a file to MinIO

        Args:
            bucket_name: Target bucket name
            object_name: Object name in MinIO (S3 key)
            file_path: Local file path to upload
            metadata: Optional metadata dictionary
            content_type: MIME type

        Returns:
            True if successful, False otherwise
        """
        try:
            self.client.fput_object(
                bucket_name,
                object_name,
                file_path,
                content_type=content_type,
                metadata=metadata
            )

            file_size = Path(file_path).stat().st_size
            logger.info(f"Uploaded {file_path} to {bucket_name}/{object_name} ({file_size} bytes)")
            return True

        except S3Error as e:
            logger.error(f"Upload failed for {file_path}: {e}")
            return False

    def upload_fastq_pair(
        self,
        user_id: str,
        job_id: str,
        sample_id: str,
        r1_path: str,
        r2_path: str
    ) -> Tuple[bool, str, str]:
        """
        Upload paired-end FASTQ files

        Args:
            user_id: User identifier
            job_id: Job identifier
            sample_id: Sample identifier
            r1_path: Path to R1 FASTQ file
            r2_path: Path to R2 FASTQ file

        Returns:
            (success, r1_object_name, r2_object_name)
        """
        bucket = self.config.BUCKET_RAW_DATA

        # Object names with hierarchical structure
        r1_object = f"{user_id}/{job_id}/{sample_id}_R1.fastq.gz"
        r2_object = f"{user_id}/{job_id}/{sample_id}_R2.fastq.gz"

        metadata = {
            "user_id": user_id,
            "job_id": job_id,
            "sample_id": sample_id,
            "upload_date": datetime.utcnow().isoformat()
        }

        r1_success = self.upload_file(bucket, r1_object, r1_path, metadata)
        r2_success = self.upload_file(bucket, r2_object, r2_path, metadata)

        return (r1_success and r2_success, r1_object, r2_object)

    def upload_results(
        self,
        user_id: str,
        job_id: str,
        sample_id: str,
        vcf_path: str,
        annotated_path: str
    ) -> bool:
        """
        Upload analysis results (VCF and annotated variants)

        Args:
            user_id: User identifier
            job_id: Job identifier
            sample_id: Sample identifier
            vcf_path: Path to VCF file
            annotated_path: Path to annotated variants file

        Returns:
            True if successful, False otherwise
        """
        bucket = self.config.BUCKET_RESULTS

        metadata = {
            "user_id": user_id,
            "job_id": job_id,
            "sample_id": sample_id,
            "analysis_date": datetime.utcnow().isoformat()
        }

        vcf_object = f"{user_id}/{job_id}/{sample_id}.vcf.gz"
        annotated_object = f"{user_id}/{job_id}/{sample_id}_annotated.txt"

        vcf_success = self.upload_file(bucket, vcf_object, vcf_path, metadata, "application/gzip")
        ann_success = self.upload_file(bucket, annotated_object, annotated_path, metadata, "text/tab-separated-values")

        return vcf_success and ann_success

    # ==================== Download Operations ====================

    def download_file(
        self,
        bucket_name: str,
        object_name: str,
        file_path: str
    ) -> bool:
        """
        Download a file from MinIO

        Args:
            bucket_name: Source bucket name
            object_name: Object name in MinIO
            file_path: Local file path to save

        Returns:
            True if successful, False otherwise
        """
        try:
            self.client.fget_object(bucket_name, object_name, file_path)
            logger.info(f"Downloaded {bucket_name}/{object_name} to {file_path}")
            return True

        except S3Error as e:
            logger.error(f"Download failed for {object_name}: {e}")
            return False

    def get_object_stream(
        self,
        bucket_name: str,
        object_name: str
    ):
        """
        Get object as a stream (for large files)

        Args:
            bucket_name: Source bucket name
            object_name: Object name in MinIO

        Returns:
            Response object (stream)
        """
        try:
            response = self.client.get_object(bucket_name, object_name)
            return response
        except S3Error as e:
            logger.error(f"Failed to get object stream for {object_name}: {e}")
            return None

    # ==================== URL Generation ====================

    def get_presigned_url(
        self,
        bucket_name: str,
        object_name: str,
        expires: timedelta = timedelta(hours=1)
    ) -> Optional[str]:
        """
        Generate presigned URL for temporary access

        Args:
            bucket_name: Bucket name
            object_name: Object name
            expires: URL expiration time (default: 1 hour)

        Returns:
            Presigned URL string or None if failed
        """
        try:
            url = self.client.presigned_get_object(
                bucket_name,
                object_name,
                expires=expires
            )
            logger.info(f"Generated presigned URL for {bucket_name}/{object_name}")
            return url

        except S3Error as e:
            logger.error(f"Failed to generate presigned URL: {e}")
            return None

    def get_results_download_urls(
        self,
        user_id: str,
        job_id: str,
        sample_id: str,
        expires_hours: int = 24
    ) -> Dict[str, str]:
        """
        Get download URLs for all result files

        Args:
            user_id: User identifier
            job_id: Job identifier
            sample_id: Sample identifier
            expires_hours: URL expiration in hours

        Returns:
            Dictionary of file_type -> presigned_url
        """
        bucket = self.config.BUCKET_RESULTS
        expires = timedelta(hours=expires_hours)

        urls = {}

        # VCF file
        vcf_object = f"{user_id}/{job_id}/{sample_id}.vcf.gz"
        urls["vcf"] = self.get_presigned_url(bucket, vcf_object, expires)

        # Annotated variants
        annotated_object = f"{user_id}/{job_id}/{sample_id}_annotated.txt"
        urls["annotated"] = self.get_presigned_url(bucket, annotated_object, expires)

        return urls

    # ==================== List Operations ====================

    def list_objects(
        self,
        bucket_name: str,
        prefix: str = "",
        recursive: bool = True
    ) -> List[Dict]:
        """
        List objects in a bucket

        Args:
            bucket_name: Bucket name
            prefix: Object prefix filter
            recursive: Recursive listing

        Returns:
            List of object dictionaries
        """
        try:
            objects = self.client.list_objects(
                bucket_name,
                prefix=prefix,
                recursive=recursive
            )

            result = []
            for obj in objects:
                result.append({
                    "name": obj.object_name,
                    "size": obj.size,
                    "last_modified": obj.last_modified,
                    "etag": obj.etag
                })

            return result

        except S3Error as e:
            logger.error(f"Failed to list objects in {bucket_name}: {e}")
            return []

    def list_user_jobs(self, user_id: str) -> List[str]:
        """
        List all jobs for a user

        Args:
            user_id: User identifier

        Returns:
            List of job IDs
        """
        bucket = self.config.BUCKET_RESULTS
        objects = self.list_objects(bucket, prefix=f"{user_id}/", recursive=False)

        # Extract unique job IDs from object paths
        job_ids = set()
        for obj in objects:
            parts = obj["name"].split("/")
            if len(parts) >= 2:
                job_ids.add(parts[1])

        return sorted(list(job_ids))

    # ==================== Delete Operations ====================

    def delete_object(
        self,
        bucket_name: str,
        object_name: str
    ) -> bool:
        """
        Delete an object from MinIO

        Args:
            bucket_name: Bucket name
            object_name: Object name

        Returns:
            True if successful, False otherwise
        """
        try:
            self.client.remove_object(bucket_name, object_name)
            logger.info(f"Deleted {bucket_name}/{object_name}")
            return True

        except S3Error as e:
            logger.error(f"Failed to delete {object_name}: {e}")
            return False

    def delete_job_data(self, user_id: str, job_id: str) -> bool:
        """
        Delete all data for a job across all buckets

        Args:
            user_id: User identifier
            job_id: Job identifier

        Returns:
            True if successful, False otherwise
        """
        buckets = [
            self.config.BUCKET_RAW_DATA,
            self.config.BUCKET_INTERMEDIATE,
            self.config.BUCKET_RESULTS,
            self.config.BUCKET_LOGS
        ]

        success = True
        prefix = f"{user_id}/{job_id}/"

        for bucket in buckets:
            try:
                objects = self.list_objects(bucket, prefix=prefix)
                for obj in objects:
                    if not self.delete_object(bucket, obj["name"]):
                        success = False
            except Exception as e:
                logger.error(f"Failed to delete job data from {bucket}: {e}")
                success = False

        return success

    # ==================== Utility Operations ====================

    def get_object_metadata(
        self,
        bucket_name: str,
        object_name: str
    ) -> Optional[Dict]:
        """
        Get object metadata

        Args:
            bucket_name: Bucket name
            object_name: Object name

        Returns:
            Metadata dictionary or None
        """
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
            logger.error(f"Failed to get metadata for {object_name}: {e}")
            return None

    def object_exists(
        self,
        bucket_name: str,
        object_name: str
    ) -> bool:
        """
        Check if object exists

        Args:
            bucket_name: Bucket name
            object_name: Object name

        Returns:
            True if exists, False otherwise
        """
        try:
            self.client.stat_object(bucket_name, object_name)
            return True
        except S3Error:
            return False

    def get_bucket_usage(self, bucket_name: str) -> Dict:
        """
        Get bucket usage statistics

        Args:
            bucket_name: Bucket name

        Returns:
            Dictionary with usage stats
        """
        objects = self.list_objects(bucket_name)

        total_size = sum(obj["size"] for obj in objects)
        object_count = len(objects)

        return {
            "bucket": bucket_name,
            "object_count": object_count,
            "total_size_bytes": total_size,
            "total_size_gb": round(total_size / (1024**3), 2)
        }

    def get_all_bucket_usage(self) -> List[Dict]:
        """
        Get usage statistics for all WES buckets

        Returns:
            List of usage dictionaries
        """
        buckets = [
            self.config.BUCKET_RAW_DATA,
            self.config.BUCKET_INTERMEDIATE,
            self.config.BUCKET_RESULTS,
            self.config.BUCKET_ARCHIVES,
            self.config.BUCKET_REFERENCE,
            self.config.BUCKET_LOGS
        ]

        return [self.get_bucket_usage(bucket) for bucket in buckets]


# Convenience function for FastAPI dependency injection
def get_minio_service() -> MinIOService:
    """
    Get MinIO service instance (for FastAPI dependency injection)

    Usage in FastAPI:
        @app.get("/files")
        def list_files(minio: MinIOService = Depends(get_minio_service)):
            return minio.list_objects("wes-results")
    """
    return MinIOService()
