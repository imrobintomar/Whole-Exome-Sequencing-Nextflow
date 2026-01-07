"""
Audit Logging Service
Tracks all sensitive operations for compliance
"""

from database_extensions import AuditLog
import json
from datetime import datetime

class AuditService:
    def __init__(self, db):
        self.db = db

    def log_action(
        self,
        action: str,
        user_id: str = None,
        resource_type: str = None,
        resource_id: str = None,
        ip_address: str = None,
        user_agent: str = None,
        metadata: dict = None
    ):
        """
        Log an audit event

        Example actions:
        - job_submit
        - job_cancel
        - subscription_created
        - subscription_canceled
        - admin_login
        - chat_message_sent
        - user_upgraded
        """
        log = AuditLog(
            user_id=user_id,
            action=action,
            resource_type=resource_type,
            resource_id=resource_id,
            ip_address=ip_address,
            user_agent=user_agent,
            metadata_json=json.dumps(metadata) if metadata else None
        )

        self.db.add(log)
        self.db.commit()
        return log

    def get_user_actions(self, user_id: str, limit: int = 100):
        """Get recent actions by user"""
        return self.db.query(AuditLog).filter(
            AuditLog.user_id == user_id
        ).order_by(AuditLog.created_at.desc()).limit(limit).all()

    def get_resource_history(self, resource_type: str, resource_id: str):
        """Get all actions on a specific resource"""
        return self.db.query(AuditLog).filter(
            AuditLog.resource_type == resource_type,
            AuditLog.resource_id == resource_id
        ).order_by(AuditLog.created_at.asc()).all()
