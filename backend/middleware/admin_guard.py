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

class AdminContext:
    """Simple data class to hold admin information without SQLAlchemy session"""
    def __init__(self, admin_user):
        self.id = admin_user.id
        self.firebase_uid = admin_user.firebase_uid
        self.email = admin_user.email
        self.is_super_admin = admin_user.is_super_admin
        self.last_login = admin_user.last_login
        self.created_at = admin_user.created_at


async def require_admin(current_user = Depends(get_current_user)):
    """
    Dependency that enforces admin-only access
    Raises 403 if user is not admin

    Admin users are identified by:
    1. Entry in admin_users table, OR
    2. Firebase UID in ADMIN_USER_UIDS environment variable

    Returns: AdminContext object (not bound to any session)
    """
    db = SessionLocal()
    try:
        # Check database first
        admin = db.query(AdminUser).filter(
            AdminUser.firebase_uid == current_user.firebase_uid
        ).first()

        if admin:
            # Update last login
            admin.last_login = datetime.utcnow()
            db.commit()

            # Return a simple data object that's not bound to the session
            return AdminContext(admin)

        # Fallback: Check environment variable
        admin_uids_str = settings.ADMIN_USER_UIDS or ""
        admin_uids = [uid.strip() for uid in admin_uids_str.split(",") if uid.strip()]

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

            # Return a simple data object that's not bound to the session
            return AdminContext(admin)

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
