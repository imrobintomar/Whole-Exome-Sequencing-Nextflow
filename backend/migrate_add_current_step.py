"""
Migration script to add current_step column to jobs table
Run this script once to update the database schema
"""
import sqlite3
from pathlib import Path

# Get database path from config
try:
    from config import settings
    db_path = settings.DATABASE_URL.replace('sqlite:///', '')
except:
    # Fallback to default path
    db_path = './wes_pipeline.db'

def migrate():
    """Add current_step column to jobs table if it doesn't exist"""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    try:
        # Check if column already exists
        cursor.execute("PRAGMA table_info(jobs)")
        columns = [column[1] for column in cursor.fetchall()]

        if 'current_step' not in columns:
            print("Adding current_step column to jobs table...")
            cursor.execute("ALTER TABLE jobs ADD COLUMN current_step TEXT")
            conn.commit()
            print("✓ Migration completed successfully!")
        else:
            print("✓ current_step column already exists, no migration needed.")

    except Exception as e:
        print(f"✗ Migration failed: {e}")
        conn.rollback()
    finally:
        conn.close()

if __name__ == "__main__":
    print("Running database migration...")
    migrate()
