"""
Activity & Audit Logging Service - Track all admin actions
"""
from datetime import datetime, timezone
from typing import Dict, List, Optional
from sqlalchemy import Column, Integer, String, DateTime, Text, JSON, desc, or_
from database import Base, SessionLocal


class ActivityLog(Base):
    """Activity log records for audit trail"""
    __tablename__ = "activity_logs"

    id = Column(Integer, primary_key=True, index=True)
    admin_uid = Column(String, index=True, nullable=False)  # Admin who performed action
    admin_email = Column(String, nullable=True)
    action_type = Column(String, index=True, nullable=False)  # ban_user, update_subscription, etc.
    resource_type = Column(String, nullable=True)  # user, job, system
    resource_id = Column(String, nullable=True)  # UID or ID of affected resource
    description = Column(Text, nullable=False)
    metadata = Column(JSON, nullable=True)  # Additional context
    ip_address = Column(String, nullable=True)
    created_at = Column(DateTime, default=lambda: datetime.now(timezone.utc), index=True)


class ActivityLogger:
    """Log and retrieve admin activity for audit purposes"""

    @staticmethod
    def log_action(
        admin_uid: str,
        action_type: str,
        description: str,
        resource_type: Optional[str] = None,
        resource_id: Optional[str] = None,
        metadata: Optional[Dict] = None,
        admin_email: Optional[str] = None,
        ip_address: Optional[str] = None
    ) -> ActivityLog:
        """
        Log an admin action

        Args:
            admin_uid: Firebase UID of admin performing action
            action_type: Type of action (ban_user, update_subscription, etc.)
            description: Human-readable description
            resource_type: Type of resource affected (user, job, system)
            resource_id: ID of affected resource
            metadata: Additional context as dictionary
            admin_email: Email of admin (optional)
            ip_address: IP address of admin (optional)
        """
        db = SessionLocal()
        try:
            log = ActivityLog(
                admin_uid=admin_uid,
                admin_email=admin_email,
                action_type=action_type,
                resource_type=resource_type,
                resource_id=resource_id,
                description=description,
                metadata=metadata,
                ip_address=ip_address
            )
            db.add(log)
            db.commit()
            db.refresh(log)

            print(f"📝 Activity Log: {action_type} by {admin_email or admin_uid} - {description}")

            return log
        finally:
            db.close()

    @staticmethod
    def get_activity_logs(
        limit: int = 100,
        offset: int = 0,
        admin_uid: Optional[str] = None,
        action_type: Optional[str] = None,
        resource_type: Optional[str] = None,
        search: Optional[str] = None,
        date_from: Optional[str] = None,
        date_to: Optional[str] = None
    ) -> Dict:
        """
        Get activity logs with filtering

        Args:
            limit: Maximum number of logs to return
            offset: Offset for pagination
            admin_uid: Filter by admin UID
            action_type: Filter by action type
            resource_type: Filter by resource type
            search: Search in description and resource_id
            date_from: Filter logs from this date
            date_to: Filter logs until this date
        """
        db = SessionLocal()
        try:
            query = db.query(ActivityLog)

            # Apply filters
            if admin_uid:
                query = query.filter(ActivityLog.admin_uid == admin_uid)

            if action_type:
                query = query.filter(ActivityLog.action_type == action_type)

            if resource_type:
                query = query.filter(ActivityLog.resource_type == resource_type)

            if search:
                query = query.filter(
                    or_(
                        ActivityLog.description.ilike(f"%{search}%"),
                        ActivityLog.resource_id.ilike(f"%{search}%"),
                        ActivityLog.admin_email.ilike(f"%{search}%")
                    )
                )

            if date_from:
                query = query.filter(ActivityLog.created_at >= date_from)

            if date_to:
                query = query.filter(ActivityLog.created_at <= date_to)

            # Get total count
            total = query.count()

            # Get paginated results
            logs = query.order_by(desc(ActivityLog.created_at)).limit(limit).offset(offset).all()

            return {
                "logs": [ActivityLogger._log_to_dict(log) for log in logs],
                "total": total,
                "limit": limit,
                "offset": offset
            }
        finally:
            db.close()

    @staticmethod
    def _log_to_dict(log: ActivityLog) -> Dict:
        """Convert activity log to dictionary"""
        return {
            "id": log.id,
            "admin_uid": log.admin_uid,
            "admin_email": log.admin_email,
            "action_type": log.action_type,
            "resource_type": log.resource_type,
            "resource_id": log.resource_id,
            "description": log.description,
            "metadata": log.metadata,
            "ip_address": log.ip_address,
            "created_at": log.created_at.isoformat() if log.created_at else None
        }

    @staticmethod
    def get_action_types() -> List[str]:
        """Get all unique action types for filtering"""
        db = SessionLocal()
        try:
            action_types = db.query(ActivityLog.action_type).distinct().all()
            return [at[0] for at in action_types if at[0]]
        finally:
            db.close()

    @staticmethod
    def get_admin_activity_summary(admin_uid: str, days: int = 30) -> Dict:
        """Get activity summary for a specific admin"""
        db = SessionLocal()
        try:
            from datetime import timedelta

            since = datetime.now(timezone.utc) - timedelta(days=days)

            logs = db.query(ActivityLog).filter(
                ActivityLog.admin_uid == admin_uid,
                ActivityLog.created_at >= since
            ).all()

            action_counts = {}
            for log in logs:
                action_counts[log.action_type] = action_counts.get(log.action_type, 0) + 1

            return {
                "admin_uid": admin_uid,
                "total_actions": len(logs),
                "action_counts": action_counts,
                "period_days": days,
                "recent_logs": [ActivityLogger._log_to_dict(log) for log in logs[:10]]
            }
        finally:
            db.close()
