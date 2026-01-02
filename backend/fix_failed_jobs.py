"""
Fix failed jobs that have NULL current_step
This script updates failed jobs to show "Initializing" instead of NULL
"""
import sqlite3

# Get database path from config
try:
    from config import settings
    db_path = settings.DATABASE_URL.replace('sqlite:///', '')
except:
    # Fallback to default path
    db_path = './wes_pipeline.db'

def fix_failed_jobs():
    """Update failed jobs with NULL current_step to show 'Initializing'"""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    try:
        # Find failed jobs with NULL current_step
        cursor.execute("""
            SELECT job_id, sample_name
            FROM jobs
            WHERE status = 'FAILED' AND (current_step IS NULL OR current_step = '')
        """)
        jobs = cursor.fetchall()

        if not jobs:
            print("✓ No failed jobs with NULL current_step found.")
            return

        print(f"Found {len(jobs)} failed job(s) with NULL current_step:")
        for job_id, sample_name in jobs:
            print(f"  - {sample_name} (ID: {job_id})")

        # Update them to show "Initializing"
        cursor.execute("""
            UPDATE jobs
            SET current_step = 'Initializing'
            WHERE status = 'FAILED' AND (current_step IS NULL OR current_step = '')
        """)
        conn.commit()

        print(f"✓ Updated {cursor.rowcount} job(s) successfully!")

    except Exception as e:
        print(f"✗ Fix failed: {e}")
        conn.rollback()
    finally:
        conn.close()

if __name__ == "__main__":
    print("Fixing failed jobs with NULL current_step...")
    fix_failed_jobs()
