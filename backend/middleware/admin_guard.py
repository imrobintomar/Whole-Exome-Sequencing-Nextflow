"""
Admin Guard Middleware
Enforces admin-only access to protected routes
"""

from fastapi import Depends, HTTPException, status
from firebase_auth import get_current_user
from database import SessionLocal
from database_extensions import AdminUser
from datetime import datetime
from config import settings

async def require_admin(current_user = Depends(get_current_user)):
    """
    Dependency that enforces admin-only access
    Raises 403 if user is not admin

    Admin users are identified by:
    1. Entry in admin_users table, OR
    2. Firebase UID in ADMIN_USER_UIDS environment variable
    """
    db = SessionLocal()
    try:
        # Debug logging
        print(f"🔍 Admin check for user: {current_user.firebase_uid}")

        # Check database first
        admin = db.query(AdminUser).filter(
            AdminUser.firebase_uid == current_user.firebase_uid
        ).first()

        if admin:
            # Update last login
            admin.last_login = datetime.utcnow()
            db.commit()
            print(f"✅ Admin found in database: {current_user.email}")
            return admin

        # Fallback: Check environment variable
        admin_uids_str = settings.ADMIN_USER_UIDS or ""
        admin_uids = [uid.strip() for uid in admin_uids_str.split(",") if uid.strip()]
        print(f"📋 Admin UIDs from env: {admin_uids}")
        print(f"🔑 Current user UID: {current_user.firebase_uid}")
        print(f"✓ Match: {current_user.firebase_uid in admin_uids}")

        if current_user.firebase_uid in admin_uids:
            # Auto-create admin user entry
            admin = AdminUser(
                firebase_uid=current_user.firebase_uid,
                email=current_user.email,
                is_super_admin=False,
                last_login=datetime.utcnow()
            )
            db.add(admin)
            db.commit()
            db.refresh(admin)
            return admin

        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Admin access required"
        )

    finally:
        db.close()


async def require_super_admin(admin = Depends(require_admin)):
    """
    Dependency that enforces super admin access
    Used for sensitive operations
    """
    if not admin.is_super_admin:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Super admin access required"
        )
    return admin
