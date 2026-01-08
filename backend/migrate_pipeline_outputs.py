"""
Migrate completed jobs to use new pipeline output file paths.
Updates file paths from old pipeline structure to new simplified structure:
- Old: ${sample}_final_annotated.tsv -> New: ${sample}_Final_.tsv
- Old: ${sample}.annovar.hg38_multianno.vcf -> New: ${sample}.annovar.hg38_multianno.txt
"""
import sqlite3
import os
from pathlib import Path

# Get database path from config
try:
    from config import settings
    db_path = settings.DATABASE_URL.replace('sqlite:///', '')
    output_dir = Path(settings.RESULTS_DIR)
except:
    # Fallback to default paths
    db_path = './wes_pipeline.db'
    output_dir = Path('./results')

def migrate_completed_jobs():
    """Update completed jobs with new file paths"""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    try:
        # Find completed jobs
        cursor.execute("""
            SELECT job_id, sample_name, annotated_vcf_path, filtered_tsv_path
            FROM jobs
            WHERE status = 'COMPLETED'
        """)
        jobs = cursor.fetchall()

        if not jobs:
            print("✓ No completed jobs found.")
            return

        print(f"Found {len(jobs)} completed job(s). Checking for file updates...")

        updated_count = 0
        for job_id, sample_name, old_annotated_path, old_tsv_path in jobs:
            needs_update = False
            new_annotated_path = old_annotated_path
            new_tsv_path = old_tsv_path

            # Search for ANNOVAR TXT file
            if not old_annotated_path or old_annotated_path.endswith('.vcf') or not os.path.exists(old_annotated_path):
                # Try multiple possible locations
                search_paths = [
                    output_dir / job_id / "output" / f"{sample_name}.annovar.hg38_multianno.txt",
                    output_dir / job_id / "output" / "Germline_VCF" / f"{sample_name}.annovar.hg38_multianno.txt",
                    output_dir / job_id / "Germline_VCF" / f"{sample_name}.annovar.hg38_multianno.txt"
                ]

                for search_path in search_paths:
                    if search_path.exists():
                        new_annotated_path = str(search_path)
                        needs_update = True
                        print(f"  ✓ {sample_name}: Found ANNOVAR TXT file at {search_path}")
                        break

            # Search for Final_.tsv file
            if not old_tsv_path or 'final_annotated.tsv' in old_tsv_path or not os.path.exists(old_tsv_path):
                # Try multiple possible locations
                search_paths = [
                    output_dir / job_id / "output" / f"{sample_name}_Final_.tsv",
                    output_dir / job_id / "Germline_VCF" / f"{sample_name}_Final_.tsv"
                ]

                for search_path in search_paths:
                    if search_path.exists():
                        new_tsv_path = str(search_path)
                        needs_update = True
                        print(f"  ✓ {sample_name}: Found Final_.tsv file at {search_path}")
                        break

            # Update if needed
            if needs_update:
                cursor.execute("""
                    UPDATE jobs
                    SET annotated_vcf_path = ?, filtered_tsv_path = ?
                    WHERE job_id = ?
                """, (new_annotated_path, new_tsv_path, job_id))
                updated_count += 1

        conn.commit()

        if updated_count > 0:
            print(f"\n✓ Successfully updated {updated_count} job(s) with new file paths!")
        else:
            print("\n✓ All jobs already have correct file paths or files not found.")

    except Exception as e:
        print(f"✗ Migration failed: {e}")
        conn.rollback()
    finally:
        conn.close()

if __name__ == "__main__":
    print("Migrating completed jobs to new pipeline output structure...")
    print("=" * 60)
    migrate_completed_jobs()
