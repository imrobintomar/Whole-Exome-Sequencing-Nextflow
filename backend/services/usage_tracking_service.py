"""
Usage Tracking Service
Manages monthly job usage limits
"""

from database_extensions import UsageTracking, Subscription, SubscriptionPlan
from datetime import datetime
from sqlalchemy import func

class UsageTrackingService:
    def __init__(self, db):
        self.db = db

    def increment_job_count(self, user_id: str):
        """
        Increment job count for current month
        Thread-safe with row locking
        """
        current_month = int(datetime.now().strftime('%Y%m'))

        # Use row locking to prevent race conditions
        usage = self.db.query(UsageTracking).filter(
            UsageTracking.user_id == user_id,
            UsageTracking.month == current_month
        ).with_for_update().first()

        if usage:
            usage.jobs_executed += 1
            usage.updated_at = datetime.utcnow()
        else:
            # Create new usage record
            subscription = self.db.query(Subscription).filter(
                Subscription.user_id == user_id,
                Subscription.status.in_(['active', 'trialing'])
            ).first()

            if subscription:
                plan = self.db.query(SubscriptionPlan).filter(
                    SubscriptionPlan.id == subscription.plan_id
                ).first()
                jobs_limit = plan.monthly_jobs_limit
            else:
                # Free tier
                free_plan = self.db.query(SubscriptionPlan).filter(
                    SubscriptionPlan.name == 'Free'
                ).first()
                jobs_limit = free_plan.monthly_jobs_limit if free_plan else 2

            usage = UsageTracking(
                user_id=user_id,
                month=current_month,
                jobs_executed=1,
                jobs_limit=jobs_limit
            )
            self.db.add(usage)

        self.db.commit()
        return usage

    def decrement_job_count(self, user_id: str):
        """
        Decrement job count (for failed jobs or refunds)
        """
        current_month = int(datetime.now().strftime('%Y%m'))

        usage = self.db.query(UsageTracking).filter(
            UsageTracking.user_id == user_id,
            UsageTracking.month == current_month
        ).first()

        if usage and usage.jobs_executed > 0:
            usage.jobs_executed -= 1
            usage.updated_at = datetime.utcnow()
            self.db.commit()

        return usage

    def get_usage_stats(self, user_id: str):
        """
        Get usage statistics for user
        Returns current month and historical data
        """
        current_month = int(datetime.now().strftime('%Y%m'))

        current = self.db.query(UsageTracking).filter(
            UsageTracking.user_id == user_id,
            UsageTracking.month == current_month
        ).first()

        historical = self.db.query(UsageTracking).filter(
            UsageTracking.user_id == user_id,
            UsageTracking.month < current_month
        ).order_by(UsageTracking.month.desc()).limit(6).all()

        return {
            "current": current,
            "historical": historical
        }

    def get_all_users_usage(self, month: int = None):
        """
        Get usage for all users (admin only)
        """
        if not month:
            month = int(datetime.now().strftime('%Y%m'))

        return self.db.query(UsageTracking).filter(
            UsageTracking.month == month
        ).all()
