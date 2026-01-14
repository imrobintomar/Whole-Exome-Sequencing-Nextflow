#!/usr/bin/env python3
"""
Test script for MinIO Service Layer
Verifies all operations work correctly
Created: 2026-01-14
"""

from minio_service import MinIOService, MinIOConfig
import tempfile
import os
from pathlib import Path

def test_minio_service():
    """
    Test all MinIO service operations
    """
    print("=" * 60)
    print("MinIO Service Layer Test")
    print("=" * 60)

    # Initialize service
    print("\n1. Initializing MinIO service...")
    try:
        minio = MinIOService()
        print("   ✓ Service initialized successfully")
    except Exception as e:
        print(f"   ✗ Failed to initialize: {e}")
        return 1

    # Test 1: List buckets via usage stats
    print("\n2. Testing bucket access...")
    try:
        usage_stats = minio.get_all_bucket_usage()
        print(f"   ✓ Found {len(usage_stats)} buckets")
        for stat in usage_stats:
            print(f"     - {stat['bucket']}: {stat['object_count']} objects, {stat['total_size_gb']} GB")
    except Exception as e:
        print(f"   ✗ Failed to get bucket usage: {e}")
        return 1

    # Test 2: Upload test file
    print("\n3. Testing file upload...")
    test_user = "test-user"
    test_job = "test-job-001"
    test_sample = "TEST001"

    try:
        # Create temporary test file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tmp:
            tmp.write("This is a test file for MinIO service layer\n")
            tmp.write(f"User: {test_user}\n")
            tmp.write(f"Job: {test_job}\n")
            tmp.write(f"Sample: {test_sample}\n")
            tmp_path = tmp.name

        # Upload to results bucket
        bucket = minio.config.BUCKET_RESULTS
        object_name = f"{test_user}/{test_job}/{test_sample}_test.txt"

        success = minio.upload_file(
            bucket,
            object_name,
            tmp_path,
            metadata={"test": "true", "user": test_user}
        )

        if success:
            print(f"   ✓ Uploaded test file to {bucket}/{object_name}")
        else:
            print("   ✗ Upload failed")
            return 1

        # Clean up temp file
        os.unlink(tmp_path)

    except Exception as e:
        print(f"   ✗ Upload test failed: {e}")
        return 1

    # Test 3: Generate presigned URL
    print("\n4. Testing presigned URL generation...")
    try:
        url = minio.get_presigned_url(bucket, object_name)
        if url:
            print(f"   ✓ Generated presigned URL")
            print(f"     URL (first 80 chars): {url[:80]}...")
        else:
            print("   ✗ Failed to generate URL")
            return 1
    except Exception as e:
        print(f"   ✗ URL generation failed: {e}")
        return 1

    # Test 4: Check object exists
    print("\n5. Testing object existence check...")
    try:
        exists = minio.object_exists(bucket, object_name)
        if exists:
            print(f"   ✓ Object exists: {object_name}")
        else:
            print("   ✗ Object not found")
            return 1
    except Exception as e:
        print(f"   ✗ Existence check failed: {e}")
        return 1

    # Test 5: Get object metadata
    print("\n6. Testing metadata retrieval...")
    try:
        metadata = minio.get_object_metadata(bucket, object_name)
        if metadata:
            print(f"   ✓ Retrieved metadata:")
            print(f"     Size: {metadata['size']} bytes")
            print(f"     Content-Type: {metadata['content_type']}")
            print(f"     Custom metadata: {metadata.get('metadata', {})}")
        else:
            print("   ✗ Failed to get metadata")
            return 1
    except Exception as e:
        print(f"   ✗ Metadata retrieval failed: {e}")
        return 1

    # Test 6: List objects
    print("\n7. Testing object listing...")
    try:
        objects = minio.list_objects(bucket, prefix=f"{test_user}/{test_job}/")
        print(f"   ✓ Found {len(objects)} objects for {test_user}/{test_job}")
        for obj in objects:
            print(f"     - {obj['name']} ({obj['size']} bytes)")
    except Exception as e:
        print(f"   ✗ Object listing failed: {e}")
        return 1

    # Test 7: Download test file
    print("\n8. Testing file download...")
    try:
        with tempfile.NamedTemporaryFile(suffix='.txt', delete=False) as tmp:
            download_path = tmp.name

        success = minio.download_file(bucket, object_name, download_path)

        if success:
            # Verify content
            with open(download_path, 'r') as f:
                content = f.read()
                if test_user in content and test_job in content:
                    print(f"   ✓ Downloaded and verified file content")
                else:
                    print("   ✗ Downloaded file content mismatch")
                    return 1

            # Clean up
            os.unlink(download_path)
        else:
            print("   ✗ Download failed")
            return 1

    except Exception as e:
        print(f"   ✗ Download test failed: {e}")
        return 1

    # Test 8: Delete test file (cleanup)
    print("\n9. Testing file deletion (cleanup)...")
    try:
        success = minio.delete_object(bucket, object_name)
        if success:
            print(f"   ✓ Deleted test file")

            # Verify deletion
            exists = minio.object_exists(bucket, object_name)
            if not exists:
                print(f"   ✓ Verified file deletion")
            else:
                print("   ⚠ File still exists after deletion")
        else:
            print("   ✗ Deletion failed")
            return 1

    except Exception as e:
        print(f"   ✗ Deletion test failed: {e}")
        return 1

    # Summary
    print("\n" + "=" * 60)
    print("✓ All tests passed successfully!")
    print("=" * 60)

    print("\n📊 Final Bucket Usage:")
    usage_stats = minio.get_all_bucket_usage()
    for stat in usage_stats:
        print(f"  {stat['bucket']:20} {stat['object_count']:5} objects  {stat['total_size_gb']:8.2f} GB")

    return 0

if __name__ == "__main__":
    import sys
    sys.exit(test_minio_service())
