"""
Analytics Service - Generate insights and trends from platform data
"""
from datetime import datetime, timezone, timedelta
from typing import Dict, List
from sqlalchemy import func, extract, desc
from database import SessionLocal, User, Job, JobStatus
from modules.saas.models import Subscription, UsageTracking


class AnalyticsService:
    """Generate analytics and insights for admin dashboard"""

    @staticmethod
    def get_jobs_over_time(days: int = 30, granularity: str = "day") -> Dict:
        """
        Get job submission trends over time

        Args:
            days: Number of days to look back
            granularity: 'day', 'week', or 'month'
        """
        db = SessionLocal()
        try:
            since = datetime.now(timezone.utc) - timedelta(days=days)

            if granularity == "day":
                # Group by day
                results = db.query(
                    func.date(Job.created_at).label('date'),
                    func.count(Job.id).label('count'),
                    Job.status
                ).filter(
                    Job.created_at >= since
                ).group_by(
                    func.date(Job.created_at),
                    Job.status
                ).order_by(
                    func.date(Job.created_at)
                ).all()

                # Organize by date and status
                data = {}
                for date, count, status in results:
                    date_str = date.isoformat() if hasattr(date, 'isoformat') else str(date)
                    if date_str not in data:
                        data[date_str] = {"total": 0, "pending": 0, "running": 0, "completed": 0, "failed": 0}
                    data[date_str][status.value if hasattr(status, 'value') else status] = count
                    data[date_str]["total"] += count

                return {
                    "granularity": "day",
                    "period_days": days,
                    "data": [{"date": k, **v} for k, v in sorted(data.items())]
                }
            else:
                # For week/month, similar logic
                return {"granularity": granularity, "data": []}

        finally:
            db.close()

    @staticmethod
    def get_user_growth(days: int = 90) -> Dict:
        """Get user registration trends"""
        db = SessionLocal()
        try:
            since = datetime.now(timezone.utc) - timedelta(days=days)

            # Group by day
            results = db.query(
                func.date(User.created_at).label('date'),
                func.count(User.id).label('count')
            ).filter(
                User.created_at >= since
            ).group_by(
                func.date(User.created_at)
            ).order_by(
                func.date(User.created_at)
            ).all()

            data = []
            cumulative = db.query(func.count(User.id)).filter(User.created_at < since).scalar() or 0

            for date, count in results:
                cumulative += count
                data.append({
                    "date": date.isoformat() if hasattr(date, 'isoformat') else str(date),
                    "new_users": count,
                    "total_users": cumulative
                })

            return {
                "period_days": days,
                "data": data,
                "total_growth": data[-1]["total_users"] - data[0]["total_users"] if data else 0
            }
        finally:
            db.close()

    @staticmethod
    def get_revenue_analytics(months: int = 12) -> Dict:
        """Get revenue trends and MRR growth"""
        db = SessionLocal()
        try:
            from modules.saas.models import SubscriptionPlan

            # Get current MRR
            active_subs = db.query(Subscription).filter(
                Subscription.status == 'active'
            ).all()

            current_mrr = 0
            plan_breakdown = {}

            for sub in active_subs:
                plan = db.query(SubscriptionPlan).filter(
                    SubscriptionPlan.id == sub.plan_id
                ).first()
                if plan:
                    current_mrr += plan.price_usd
                    plan_breakdown[plan.name] = plan_breakdown.get(plan.name, 0) + 1

            # Get subscription growth over time
            since = datetime.now(timezone.utc) - timedelta(days=months * 30)

            monthly_revenue = db.query(
                extract('year', Subscription.created_at).label('year'),
                extract('month', Subscription.created_at).label('month'),
                func.count(Subscription.id).label('new_subscriptions')
            ).filter(
                Subscription.created_at >= since
            ).group_by(
                extract('year', Subscription.created_at),
                extract('month', Subscription.created_at)
            ).order_by(
                extract('year', Subscription.created_at),
                extract('month', Subscription.created_at)
            ).all()

            revenue_data = []
            for year, month, new_subs in monthly_revenue:
                revenue_data.append({
                    "year": int(year),
                    "month": int(month),
                    "new_subscriptions": new_subs
                })

            return {
                "current_mrr": current_mrr,
                "active_subscriptions": len(active_subs),
                "plan_breakdown": plan_breakdown,
                "monthly_data": revenue_data
            }
        finally:
            db.close()

    @staticmethod
    def get_job_success_rates(days: int = 30) -> Dict:
        """Calculate job success and failure rates"""
        db = SessionLocal()
        try:
            since = datetime.now(timezone.utc) - timedelta(days=days)

            # Get job counts by status
            status_counts = db.query(
                Job.status,
                func.count(Job.id).label('count')
            ).filter(
                Job.created_at >= since
            ).group_by(
                Job.status
            ).all()

            total = sum(count for _, count in status_counts)
            status_data = {}

            for status, count in status_counts:
                status_key = status.value if hasattr(status, 'value') else status
                status_data[status_key] = {
                    "count": count,
                    "percentage": round((count / total * 100), 2) if total > 0 else 0
                }

            # Calculate success rate
            completed = status_data.get("completed", {}).get("count", 0)
            failed = status_data.get("failed", {}).get("count", 0)

            success_rate = round((completed / (completed + failed) * 100), 2) if (completed + failed) > 0 else 0

            return {
                "period_days": days,
                "total_jobs": total,
                "status_breakdown": status_data,
                "success_rate": success_rate,
                "failure_rate": round(100 - success_rate, 2)
            }
        finally:
            db.close()

    @staticmethod
    def get_top_users_by_usage(limit: int = 10) -> Dict:
        """Get top users by job count"""
        db = SessionLocal()
        try:
            # Get job counts per user
            user_job_counts = db.query(
                Job.user_id,
                func.count(Job.id).label('job_count')
            ).group_by(
                Job.user_id
            ).order_by(
                desc('job_count')
            ).limit(limit).all()

            top_users = []
            for user_id, job_count in user_job_counts:
                user = db.query(User).filter(User.firebase_uid == user_id).first()
                if user:
                    top_users.append({
                        "user_id": user_id,
                        "email": user.email,
                        "job_count": job_count
                    })

            return {
                "top_users": top_users,
                "limit": limit
            }
        finally:
            db.close()

    @staticmethod
    def get_user_retention(cohort_days: int = 30) -> Dict:
        """Calculate user retention metrics"""
        db = SessionLocal()
        try:
            since = datetime.now(timezone.utc) - timedelta(days=cohort_days)

            # Get users created in period
            cohort_users = db.query(User).filter(
                User.created_at >= since
            ).all()

            # Check how many have submitted jobs
            active_users = 0
            for user in cohort_users:
                job_count = db.query(func.count(Job.id)).filter(
                    Job.user_id == user.firebase_uid
                ).scalar()
                if job_count > 0:
                    active_users += 1

            retention_rate = round((active_users / len(cohort_users) * 100), 2) if cohort_users else 0

            return {
                "cohort_period_days": cohort_days,
                "total_users": len(cohort_users),
                "active_users": active_users,
                "retention_rate": retention_rate
            }
        finally:
            db.close()

    @staticmethod
    def get_platform_overview() -> Dict:
        """Get comprehensive platform overview"""
        db = SessionLocal()
        try:
            total_users = db.query(func.count(User.id)).scalar() or 0
            total_jobs = db.query(func.count(Job.id)).scalar() or 0

            # Last 7 days activity
            week_ago = datetime.now(timezone.utc) - timedelta(days=7)
            new_users_week = db.query(func.count(User.id)).filter(User.created_at >= week_ago).scalar() or 0
            new_jobs_week = db.query(func.count(Job.id)).filter(Job.created_at >= week_ago).scalar() or 0

            # Active subscriptions
            active_subs = db.query(func.count(Subscription.id)).filter(
                Subscription.status == 'active'
            ).scalar() or 0

            return {
                "total_users": total_users,
                "total_jobs": total_jobs,
                "new_users_this_week": new_users_week,
                "new_jobs_this_week": new_jobs_week,
                "active_subscriptions": active_subs
            }
        finally:
            db.close()
