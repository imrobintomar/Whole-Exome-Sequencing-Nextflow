#!/usr/bin/env python3
"""
Test script for the user details endpoint
"""
import sys
from pathlib import Path

# Add backend to path
sys.path.insert(0, str(Path(__file__).parent / 'backend'))

from database import SessionLocal, User
from modules.admin.routes import get_user_details
from middleware.admin_guard import AdminContext

def test_endpoint():
    """Test the user details endpoint directly"""
    db = SessionLocal()

    try:
        # Get first user from database
        user = db.query(User).first()

        if not user:
            print("❌ No users found in database")
            return

        print(f"✅ Found user: {user.email} (UID: {user.firebase_uid})")

        # Create a mock admin context
        class MockAdmin:
            firebase_uid = "test_admin"

        # Call the endpoint
        result = get_user_details(user.firebase_uid, MockAdmin())

        print(f"\n✅ Endpoint returned successfully!")
        print(f"   User: {result['user']['email']}")
        print(f"   Jobs: {len(result['jobs'])}")
        print(f"   Notes: {len(result['notes'])}")
        print(f"   Tags: {len(result['tags'])}")
        print(f"   Payment History: {len(result['payment_history'])}")
        print(f"   Support Tickets: {len(result['support_tickets'])}")
        print(f"   Activity Timeline: {len(result['activity_timeline'])}")

    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
    finally:
        db.close()

if __name__ == "__main__":
    test_endpoint()
