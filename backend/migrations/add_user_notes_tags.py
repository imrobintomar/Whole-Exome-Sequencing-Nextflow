"""
Migration script to add user_notes and user_tags tables
Run this script to create the new tables for the user details feature
"""

import sys
from pathlib import Path

# Add parent directory to path to import database modules
sys.path.insert(0, str(Path(__file__).parent.parent))

from database import Base, engine
from database_extensions import UserNote, UserTag

def migrate():
    """Create user_notes and user_tags tables"""
    print("Creating user_notes and user_tags tables...")

    # Create tables
    Base.metadata.create_all(bind=engine, tables=[UserNote.__table__, UserTag.__table__])

    print("✅ Migration completed successfully!")
    print("Tables created:")
    print("  - user_notes")
    print("  - user_tags")

if __name__ == "__main__":
    migrate()
