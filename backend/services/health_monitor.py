"""
Health Monitoring Service - Tracks system health and generates alerts
"""
import psutil
import os
from datetime import datetime, timezone, timedelta
from typing import Dict, List, Optional
from sqlalchemy import Column, Integer, String, DateTime, Boolean, Float, Text, desc
from database import Base, SessionLocal

class HealthAlert(Base):
    """Health alert records"""
    __tablename__ = "health_alerts"

    id = Column(Integer, primary_key=True, index=True)
    alert_type = Column(String, nullable=False)  # cpu, memory, disk, process
    severity = Column(String, nullable=False)  # warning, critical
    message = Column(Text, nullable=False)
    value = Column(Float, nullable=False)
    threshold = Column(Float, nullable=False)
    created_at = Column(DateTime, default=lambda: datetime.now(timezone.utc))
    resolved = Column(Boolean, default=False)
    resolved_at = Column(DateTime, nullable=True)


class HealthMonitor:
    """Monitor system health and generate alerts"""

    # Configurable thresholds - can be overridden by environment variables
    THRESHOLDS = {
        "cpu_warning": 70,
        "cpu_critical": int(os.getenv("HEALTH_CPU_THRESHOLD", "90")),
        "memory_warning": 70,
        "memory_critical": int(os.getenv("HEALTH_MEMORY_THRESHOLD", "90")),
        "disk_warning": 80,
        "disk_critical": int(os.getenv("HEALTH_DISK_THRESHOLD", "90")),
    }

    # Email notification settings
    EMAIL_ALERTS_ENABLED = os.getenv("HEALTH_EMAIL_ALERTS_ENABLED", "true").lower() == "true"

    @staticmethod
    def check_system_health() -> Dict:
        """
        Check system health and generate alerts if thresholds exceeded
        Returns current health status and any active alerts
        """
        db = SessionLocal()
        try:
            alerts = []

            # Get current metrics
            cpu = psutil.cpu_percent(interval=1)
            memory = psutil.virtual_memory().percent
            disk = psutil.disk_usage('/').percent

            # Prepare metrics for email alerts
            current_metrics = {
                "cpu": cpu,
                "memory": memory,
                "disk": disk
            }

            # Check CPU
            if cpu >= HealthMonitor.THRESHOLDS["cpu_critical"]:
                alert = HealthMonitor._create_alert(
                    db, "cpu", "critical",
                    f"CPU usage critically high at {cpu:.1f}%",
                    cpu, HealthMonitor.THRESHOLDS["cpu_critical"],
                    current_metrics
                )
                alerts.append(alert)
            elif cpu >= HealthMonitor.THRESHOLDS["cpu_warning"]:
                alert = HealthMonitor._create_alert(
                    db, "cpu", "warning",
                    f"CPU usage high at {cpu:.1f}%",
                    cpu, HealthMonitor.THRESHOLDS["cpu_warning"],
                    current_metrics
                )
                alerts.append(alert)

            # Check Memory
            if memory >= HealthMonitor.THRESHOLDS["memory_critical"]:
                alert = HealthMonitor._create_alert(
                    db, "memory", "critical",
                    f"Memory usage critically high at {memory:.1f}%",
                    memory, HealthMonitor.THRESHOLDS["memory_critical"],
                    current_metrics
                )
                alerts.append(alert)
            elif memory >= HealthMonitor.THRESHOLDS["memory_warning"]:
                alert = HealthMonitor._create_alert(
                    db, "memory", "warning",
                    f"Memory usage high at {memory:.1f}%",
                    memory, HealthMonitor.THRESHOLDS["memory_warning"],
                    current_metrics
                )
                alerts.append(alert)

            # Check Disk
            if disk >= HealthMonitor.THRESHOLDS["disk_critical"]:
                alert = HealthMonitor._create_alert(
                    db, "disk", "critical",
                    f"Disk usage critically high at {disk:.1f}%",
                    disk, HealthMonitor.THRESHOLDS["disk_critical"],
                    current_metrics
                )
                alerts.append(alert)
            elif disk >= HealthMonitor.THRESHOLDS["disk_warning"]:
                alert = HealthMonitor._create_alert(
                    db, "disk", "warning",
                    f"Disk usage high at {disk:.1f}%",
                    disk, HealthMonitor.THRESHOLDS["disk_warning"],
                    current_metrics
                )
                alerts.append(alert)

            # Get active unresolved alerts from last 24 hours
            active_alerts = db.query(HealthAlert).filter(
                HealthAlert.resolved == False,
                HealthAlert.created_at >= datetime.now(timezone.utc) - timedelta(hours=24)
            ).all()

            return {
                "status": "healthy" if len(alerts) == 0 else "unhealthy",
                "timestamp": datetime.now(timezone.utc).isoformat(),
                "metrics": {
                    "cpu": cpu,
                    "memory": memory,
                    "disk": disk
                },
                "new_alerts": [HealthMonitor._alert_to_dict(a) for a in alerts],
                "active_alerts_count": len(active_alerts),
                "active_alerts": [HealthMonitor._alert_to_dict(a) for a in active_alerts]
            }
        finally:
            db.close()

    @staticmethod
    def _create_alert(db, alert_type: str, severity: str, message: str, value: float, threshold: float, metrics: Optional[Dict] = None) -> HealthAlert:
        """Create a new health alert if similar alert doesn't exist in last 30 minutes"""
        # Check if similar alert exists in last 30 minutes (avoid spam)
        recent_alert = db.query(HealthAlert).filter(
            HealthAlert.alert_type == alert_type,
            HealthAlert.severity == severity,
            HealthAlert.created_at >= datetime.now(timezone.utc) - timedelta(minutes=30)
        ).first()

        if recent_alert:
            return recent_alert

        # Create new alert
        alert = HealthAlert(
            alert_type=alert_type,
            severity=severity,
            message=message,
            value=value,
            threshold=threshold
        )
        db.add(alert)
        db.commit()
        db.refresh(alert)

        print(f"🚨 Health Alert: {severity.upper()} - {message}")

        # Send email notification if enabled
        if HealthMonitor.EMAIL_ALERTS_ENABLED:
            try:
                from services.email_service import get_email_service
                email_service = get_email_service()
                email_service.send_health_alert(
                    alert_type=alert_type,
                    severity=severity,
                    message=message,
                    value=value,
                    threshold=threshold,
                    metrics=metrics
                )
            except Exception as e:
                print(f"⚠️ Failed to send email alert: {str(e)}")

        return alert

    @staticmethod
    def _alert_to_dict(alert: HealthAlert) -> Dict:
        """Convert alert to dictionary"""
        return {
            "id": alert.id,
            "type": alert.alert_type,
            "severity": alert.severity,
            "message": alert.message,
            "value": alert.value,
            "threshold": alert.threshold,
            "created_at": alert.created_at.isoformat() if alert.created_at else None,
            "resolved": alert.resolved,
            "resolved_at": alert.resolved_at.isoformat() if alert.resolved_at else None
        }

    @staticmethod
    def get_alert_history(limit: int = 50, severity: Optional[str] = None) -> List[Dict]:
        """Get health alert history"""
        db = SessionLocal()
        try:
            query = db.query(HealthAlert).order_by(desc(HealthAlert.created_at))

            if severity:
                query = query.filter(HealthAlert.severity == severity)

            alerts = query.limit(limit).all()
            return [HealthMonitor._alert_to_dict(a) for a in alerts]
        finally:
            db.close()

    @staticmethod
    def resolve_alert(alert_id: int) -> bool:
        """Mark an alert as resolved"""
        db = SessionLocal()
        try:
            alert = db.query(HealthAlert).filter(HealthAlert.id == alert_id).first()
            if alert:
                alert.resolved = True
                alert.resolved_at = datetime.now(timezone.utc)
                db.commit()
                return True
            return False
        finally:
            db.close()

    @staticmethod
    def update_thresholds(thresholds: Dict[str, float]) -> Dict:
        """Update alert thresholds"""
        for key, value in thresholds.items():
            if key in HealthMonitor.THRESHOLDS:
                HealthMonitor.THRESHOLDS[key] = value

        return HealthMonitor.THRESHOLDS
