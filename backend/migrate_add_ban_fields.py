#!/usr/bin/env python3
"""
Database migration: Add ban/suspend fields to users table
"""

import sqlite3
from pathlib import Path

DB_PATH = "./wes_pipeline.db"

def migrate():
    print("🔄 Running migration: Add ban/suspend fields to users table")

    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    try:
        # Check if columns already exist
        cursor.execute("PRAGMA table_info(users)")
        columns = [col[1] for col in cursor.fetchall()]

        migrations_run = []

        if 'is_active' not in columns:
            print("  Adding column: is_active")
            cursor.execute("ALTER TABLE users ADD COLUMN is_active BOOLEAN DEFAULT 1")
            migrations_run.append('is_active')

        if 'is_banned' not in columns:
            print("  Adding column: is_banned")
            cursor.execute("ALTER TABLE users ADD COLUMN is_banned BOOLEAN DEFAULT 0")
            migrations_run.append('is_banned')

        if 'ban_reason' not in columns:
            print("  Adding column: ban_reason")
            cursor.execute("ALTER TABLE users ADD COLUMN ban_reason TEXT")
            migrations_run.append('ban_reason')

        if 'banned_at' not in columns:
            print("  Adding column: banned_at")
            cursor.execute("ALTER TABLE users ADD COLUMN banned_at DATETIME")
            migrations_run.append('banned_at')

        if 'banned_by' not in columns:
            print("  Adding column: banned_by")
            cursor.execute("ALTER TABLE users ADD COLUMN banned_by TEXT")
            migrations_run.append('banned_by')

        if migrations_run:
            conn.commit()
            print(f"✅ Migration complete! Added columns: {', '.join(migrations_run)}")
        else:
            print("✅ All columns already exist. No migration needed.")

    except Exception as e:
        conn.rollback()
        print(f"❌ Migration failed: {e}")
        raise
    finally:
        conn.close()

if __name__ == "__main__":
    migrate()
