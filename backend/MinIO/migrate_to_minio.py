#!/usr/bin/env python3
"""
Migrate existing WES pipeline data to MinIO
Comprehensive data migration script
Created: 2026-01-14
"""

from minio_service import MinIOService
from pathlib import Path
import sys
import os
import time
from datetime import datetime

class DataMigration:
    """
    Handles migration of WES pipeline data to MinIO
    """

    def __init__(self):
        self.minio = MinIOService()
        self.stats = {
            "reference": {"count": 0, "size": 0, "failed": []},
            "raw_data": {"count": 0, "size": 0, "failed": []},
            "intermediate": {"count": 0, "size": 0, "failed": []},
            "results": {"count": 0, "size": 0, "failed": []},
        }
        self.start_time = time.time()

    def format_size(self, bytes_size):
        """Convert bytes to human readable format"""
        for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
            if bytes_size < 1024.0:
                return f"{bytes_size:.2f} {unit}"
            bytes_size /= 1024.0
        return f"{bytes_size:.2f} PB"

    def upload_with_progress(self, file_path, bucket, object_name, metadata=None):
        """Upload file with progress tracking"""
        file_size = Path(file_path).stat().st_size
        print(f"  Uploading: {Path(file_path).name} ({self.format_size(file_size)})")

        success = self.minio.upload_file(
            bucket,
            object_name,
            str(file_path),
            metadata=metadata
        )

        return success, file_size

    def migrate_reference_genomes(self):
        """
        Migrate reference genomes and indices to wes-reference bucket
        """
        print("\n" + "=" * 70)
        print("1. Migrating Reference Genomes")
        print("=" * 70)

        bucket = self.minio.config.BUCKET_REFERENCE
        reference_files = [
            # hg38 reference genome
            ("/media/drprabudh/m1/hg38/hg38.fa", "hg38/hg38.fa"),
            ("/media/drprabudh/m1/hg38/hg38.fa.fai", "hg38/hg38.fa.fai"),
            ("/media/drprabudh/m1/hg38/hg38.dict", "hg38/hg38.dict"),

            # BWA indices
            ("/media/drprabudh/m1/hg38/hg38.fa.amb", "hg38/hg38.fa.amb"),
            ("/media/drprabudh/m1/hg38/hg38.fa.ann", "hg38/hg38.fa.ann"),
            ("/media/drprabudh/m1/hg38/hg38.fa.bwt", "hg38/hg38.fa.bwt"),
            ("/media/drprabudh/m1/hg38/hg38.fa.pac", "hg38/hg38.fa.pac"),
            ("/media/drprabudh/m1/hg38/hg38.fa.sa", "hg38/hg38.fa.sa"),

            # Exome intervals
            ("/media/drprabudh/m1/Downloads/WES_Final_Agilent_V8_Hg38_1.bed", "intervals/WES_Final_Agilent_V8_Hg38_1.bed"),

            # Known sites VCF files
            ("/media/drprabudh/m1/vcf_file/Homo_sapiens_assembly38.known_indels.vcf.gz", "known_sites/Homo_sapiens_assembly38.known_indels.vcf.gz"),
            ("/media/drprabudh/m1/vcf_file/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi", "known_sites/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"),
            ("/media/drprabudh/m1/vcf_file/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz", "known_sites/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"),
            ("/media/drprabudh/m1/vcf_file/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi", "known_sites/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"),
        ]

        metadata = {
            "category": "reference",
            "upload_date": datetime.utcnow().isoformat()
        }

        for file_path, object_name in reference_files:
            if not Path(file_path).exists():
                print(f"  ⚠ SKIP: {file_path} (not found)")
                continue

            success, size = self.upload_with_progress(file_path, bucket, object_name, metadata)

            if success:
                self.stats["reference"]["count"] += 1
                self.stats["reference"]["size"] += size
                print(f"    ✓ Uploaded to {bucket}/{object_name}")
            else:
                self.stats["reference"]["failed"].append(file_path)
                print(f"    ✗ Failed: {file_path}")

    def migrate_annovar_databases(self):
        """
        Migrate ANNOVAR databases to wes-reference bucket
        """
        print("\n" + "=" * 70)
        print("2. Migrating ANNOVAR Databases")
        print("=" * 70)

        bucket = self.minio.config.BUCKET_REFERENCE
        annovar_base = Path("/media/drprabudh/m1/annovar")

        if not annovar_base.exists():
            print("  ⚠ ANNOVAR directory not found")
            return

        # Find all ANNOVAR database files
        annovar_files = []

        # hg38_humandb directory
        humandb_dir = annovar_base / "hg38_humandb"
        if humandb_dir.exists():
            for file in humandb_dir.rglob("*"):
                if file.is_file():
                    rel_path = file.relative_to(annovar_base)
                    annovar_files.append((str(file), f"annovar/{rel_path}"))

        # dbNSFP files
        dbnsfp_dir = annovar_base / "dbnsfp"
        if dbnsfp_dir.exists():
            for file in dbnsfp_dir.rglob("*"):
                if file.is_file():
                    rel_path = file.relative_to(annovar_base)
                    annovar_files.append((str(file), f"annovar/{rel_path}"))

        # 1000 genomes
        tg_dir = annovar_base / "1000"
        if tg_dir.exists():
            for file in tg_dir.rglob("*"):
                if file.is_file():
                    rel_path = file.relative_to(annovar_base)
                    annovar_files.append((str(file), f"annovar/{rel_path}"))

        # Perl scripts
        for script in ["table_annovar.pl", "annotate_variation.pl"]:
            script_path = annovar_base / script
            if script_path.exists():
                annovar_files.append((str(script_path), f"annovar/{script}"))

        print(f"  Found {len(annovar_files)} ANNOVAR files")

        metadata = {
            "category": "annovar",
            "upload_date": datetime.utcnow().isoformat()
        }

        for i, (file_path, object_name) in enumerate(annovar_files, 1):
            print(f"\n  [{i}/{len(annovar_files)}]", end=" ")
            success, size = self.upload_with_progress(file_path, bucket, object_name, metadata)

            if success:
                self.stats["reference"]["count"] += 1
                self.stats["reference"]["size"] += size
            else:
                self.stats["reference"]["failed"].append(file_path)
                print(f"    ✗ Failed")

    def migrate_raw_fastq_files(self):
        """
        Migrate raw FASTQ files to wes-raw-data bucket
        """
        print("\n" + "=" * 70)
        print("3. Migrating Raw FASTQ Files")
        print("=" * 70)

        bucket = self.minio.config.BUCKET_RAW_DATA
        uploads_dir = Path("/media/drprabudh/m3/Nextflow-Script/WholeExome/backend/uploads")

        if not uploads_dir.exists():
            print("  ⚠ Uploads directory not found")
            return

        # Process each job directory
        for job_dir in uploads_dir.iterdir():
            if not job_dir.is_dir():
                continue

            job_id = job_dir.name
            print(f"\n  Job: {job_id}")

            # Find FASTQ files
            fastq_files = list(job_dir.glob("*.fastq.gz")) + list(job_dir.glob("*.fq.gz"))

            metadata = {
                "job_id": job_id,
                "user_id": "migrated",  # Set actual user_id if available
                "upload_date": datetime.utcnow().isoformat()
            }

            for fastq_file in fastq_files:
                object_name = f"migrated/{job_id}/{fastq_file.name}"
                success, size = self.upload_with_progress(str(fastq_file), bucket, object_name, metadata)

                if success:
                    self.stats["raw_data"]["count"] += 1
                    self.stats["raw_data"]["size"] += size
                    print(f"    ✓ {fastq_file.name}")
                else:
                    self.stats["raw_data"]["failed"].append(str(fastq_file))
                    print(f"    ✗ Failed: {fastq_file.name}")

    def migrate_intermediate_files(self):
        """
        Migrate intermediate BAM files to wes-intermediate bucket
        """
        print("\n" + "=" * 70)
        print("4. Migrating Intermediate Files (BAM)")
        print("=" * 70)

        bucket = self.minio.config.BUCKET_INTERMEDIATE
        work_dir = Path("/media/drprabudh/m3/Nextflow-Script/WholeExome/backend/work")

        if not work_dir.exists():
            print("  ⚠ Work directory not found")
            return

        # Find all BAM files (these are large, selective migration)
        bam_files = []
        for bam_file in work_dir.rglob("*.bam"):
            # Only migrate specific types to save space
            if any(x in bam_file.name for x in ["sorted", "markdup", "recall"]):
                bam_files.append(bam_file)

        print(f"  Found {len(bam_files)} intermediate BAM files")
        print(f"  ⚠ Note: Only migrating sorted/markdup/recall BAMs to save space")

        # Ask for confirmation
        if len(bam_files) > 0:
            total_size = sum(f.stat().st_size for f in bam_files)
            print(f"  Total size: {self.format_size(total_size)}")

            response = input("\n  Migrate intermediate files? (y/N): ").strip().lower()
            if response != 'y':
                print("  Skipping intermediate file migration")
                return

        metadata = {
            "category": "intermediate",
            "upload_date": datetime.utcnow().isoformat()
        }

        for i, bam_file in enumerate(bam_files, 1):
            # Extract sample name from filename
            sample_name = bam_file.stem.split('.')[0].split('_')[0]
            object_name = f"migrated/{sample_name}/{bam_file.name}"

            print(f"\n  [{i}/{len(bam_files)}]", end=" ")
            success, size = self.upload_with_progress(str(bam_file), bucket, object_name, metadata)

            if success:
                self.stats["intermediate"]["count"] += 1
                self.stats["intermediate"]["size"] += size
            else:
                self.stats["intermediate"]["failed"].append(str(bam_file))

    def migrate_final_results(self):
        """
        Migrate final results (VCF, annotated variants) to wes-results bucket
        """
        print("\n" + "=" * 70)
        print("5. Migrating Final Results")
        print("=" * 70)

        bucket = self.minio.config.BUCKET_RESULTS
        results_dir = Path("/media/drprabudh/m3/Nextflow-Script/WholeExome/backend/results")

        if not results_dir.exists():
            print("  ⚠ Results directory not found")
            return

        # Process each job directory
        for job_dir in results_dir.iterdir():
            if not job_dir.is_dir():
                continue

            job_id = job_dir.name
            print(f"\n  Job: {job_id}")

            # Find result files
            result_files = []
            result_files.extend(job_dir.rglob("*.vcf.gz"))
            result_files.extend(job_dir.rglob("*.vcf"))
            result_files.extend(job_dir.rglob("*_Final_*.txt"))
            result_files.extend(job_dir.rglob("*multianno.txt"))

            metadata = {
                "job_id": job_id,
                "user_id": "migrated",
                "upload_date": datetime.utcnow().isoformat()
            }

            for result_file in result_files:
                object_name = f"migrated/{job_id}/{result_file.name}"
                success, size = self.upload_with_progress(str(result_file), bucket, object_name, metadata)

                if success:
                    self.stats["results"]["count"] += 1
                    self.stats["results"]["size"] += size
                    print(f"    ✓ {result_file.name}")
                else:
                    self.stats["results"]["failed"].append(str(result_file))
                    print(f"    ✗ Failed: {result_file.name}")

    def print_summary(self):
        """
        Print migration summary report
        """
        elapsed = time.time() - self.start_time

        print("\n" + "=" * 70)
        print("MIGRATION SUMMARY")
        print("=" * 70)

        total_files = 0
        total_size = 0

        categories = [
            ("Reference Genomes & DBs", "reference"),
            ("Raw FASTQ Files", "raw_data"),
            ("Intermediate Files", "intermediate"),
            ("Final Results", "results")
        ]

        for name, key in categories:
            count = self.stats[key]["count"]
            size = self.stats[key]["size"]
            failed = len(self.stats[key]["failed"])

            print(f"\n{name}:")
            print(f"  Uploaded:  {count} files ({self.format_size(size)})")
            if failed > 0:
                print(f"  Failed:    {failed} files")

            total_files += count
            total_size += size

        print("\n" + "-" * 70)
        print(f"Total:       {total_files} files ({self.format_size(total_size)})")
        print(f"Duration:    {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")
        print("=" * 70)

        # Show bucket usage
        print("\nMinIO Bucket Usage:")
        usage_stats = self.minio.get_all_bucket_usage()
        for stat in usage_stats:
            if stat['object_count'] > 0:
                print(f"  {stat['bucket']:20} {stat['object_count']:5} objects  {stat['total_size_gb']:8.2f} GB")

        # Save failed files to log
        if any(self.stats[key]["failed"] for key in self.stats):
            log_file = "/media/drprabudh/m3/Nextflow-Script/WholeExome/backend/MinIO/migration_failures.log"
            with open(log_file, 'w') as f:
                f.write(f"Migration Failures - {datetime.now().isoformat()}\n")
                f.write("=" * 70 + "\n\n")
                for key in self.stats:
                    if self.stats[key]["failed"]:
                        f.write(f"\n{key.upper()}:\n")
                        for failed_file in self.stats[key]["failed"]:
                            f.write(f"  {failed_file}\n")
            print(f"\n⚠ Failed files logged to: {log_file}")

    def run(self, skip_reference=False, skip_annovar=False, skip_intermediate=False):
        """
        Run complete migration

        Args:
            skip_reference: Skip reference genome migration
            skip_annovar: Skip ANNOVAR database migration
            skip_intermediate: Skip intermediate file migration
        """
        print("=" * 70)
        print("WES Pipeline Data Migration to MinIO")
        print("=" * 70)
        print(f"Start Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

        try:
            if not skip_reference:
                self.migrate_reference_genomes()

            if not skip_annovar:
                self.migrate_annovar_databases()

            self.migrate_raw_fastq_files()

            if not skip_intermediate:
                self.migrate_intermediate_files()

            self.migrate_final_results()

            self.print_summary()
            return 0

        except KeyboardInterrupt:
            print("\n\n⚠ Migration interrupted by user")
            self.print_summary()
            return 1

        except Exception as e:
            print(f"\n✗ Migration failed with error: {e}")
            import traceback
            traceback.print_exc()
            return 1


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Migrate WES pipeline data to MinIO")
    parser.add_argument("--skip-reference", action="store_true", help="Skip reference genome migration")
    parser.add_argument("--skip-annovar", action="store_true", help="Skip ANNOVAR database migration")
    parser.add_argument("--skip-intermediate", action="store_true", help="Skip intermediate file migration")
    parser.add_argument("--reference-only", action="store_true", help="Only migrate reference data")

    args = parser.parse_args()

    migrator = DataMigration()

    if args.reference_only:
        # Only migrate references
        migrator.migrate_reference_genomes()
        migrator.migrate_annovar_databases()
        migrator.print_summary()
    else:
        # Full migration with options
        sys.exit(migrator.run(
            skip_reference=args.skip_reference,
            skip_annovar=args.skip_annovar,
            skip_intermediate=args.skip_intermediate
        ))
