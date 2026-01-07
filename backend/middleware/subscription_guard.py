"""
Subscription Guard Middleware
Enforces subscription limits before job execution
"""

from fastapi import Depends, HTTPException, status
from firebase_auth import get_current_user
from database import SessionLocal
from database_extensions import Subscription, SubscriptionPlan, UsageTracking
from datetime import datetime

class SubscriptionService:
    def __init__(self, db):
        self.db = db

    def get_active_subscription(self, user_id: str):
        """Get user's active subscription"""
        return self.db.query(Subscription).filter(
            Subscription.user_id == user_id,
            Subscription.status.in_(['active', 'trialing'])
        ).first()

    def get_plan(self, plan_id: int):
        """Get subscription plan details"""
        return self.db.query(SubscriptionPlan).filter(
            SubscriptionPlan.id == plan_id
        ).first()

    def get_current_month_usage(self, user_id: str):
        """Get or create current month usage tracking"""
        current_month = int(datetime.now().strftime('%Y%m'))

        usage = self.db.query(UsageTracking).filter(
            UsageTracking.user_id == user_id,
            UsageTracking.month == current_month
        ).first()

        if not usage:
            # Get user's plan to set limit
            subscription = self.get_active_subscription(user_id)
            if not subscription:
                # No subscription - use free plan limits
                free_plan = self.db.query(SubscriptionPlan).filter(
                    SubscriptionPlan.name == 'Free'
                ).first()
                jobs_limit = free_plan.monthly_jobs_limit if free_plan else 2
            else:
                plan = self.get_plan(subscription.plan_id)
                jobs_limit = plan.monthly_jobs_limit

            usage = UsageTracking(
                user_id=user_id,
                month=current_month,
                jobs_executed=0,
                jobs_limit=jobs_limit
            )
            self.db.add(usage)
            self.db.commit()
            self.db.refresh(usage)

        return usage


async def enforce_usage_limit(current_user = Depends(get_current_user)):
    """
    Dependency that checks if user can submit a job
    Enforced BEFORE job submission

    Checks:
    1. User has active subscription or is on free tier
    2. User has not exceeded monthly job limit
    3. Subscription is not past_due or canceled

    Returns subscription data for logging
    """
    db = SessionLocal()
    try:
        service = SubscriptionService(db)

        # Get active subscription (may be None for free users)
        subscription = service.get_active_subscription(current_user.uid)

        # Determine plan and limits
        if subscription:
            # Paid subscription
            if subscription.status not in ['active', 'trialing']:
                raise HTTPException(
                    status_code=status.HTTP_402_PAYMENT_REQUIRED,
                    detail=f"Subscription {subscription.status}. Please update payment method."
                )
            plan = service.get_plan(subscription.plan_id)
        else:
            # Free tier - get free plan
            plan = db.query(SubscriptionPlan).filter(
                SubscriptionPlan.name == 'Free'
            ).first()

            if not plan:
                raise HTTPException(
                    status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                    detail="Free plan not configured. Please contact support."
                )

        # Check usage limit
        usage = service.get_current_month_usage(current_user.uid)

        if usage.jobs_executed >= plan.monthly_jobs_limit:
            raise HTTPException(
                status_code=status.HTTP_429_TOO_MANY_REQUESTS,
                detail=f"Monthly job limit reached ({plan.monthly_jobs_limit} jobs). " +
                       f"Upgrade your plan or wait for next billing cycle."
            )

        return {
            "subscription": subscription,
            "usage": usage,
            "plan": plan
        }
    finally:
        db.close()


async def check_chat_access(current_user = Depends(get_current_user)):
    """
    Check if user has access to chat support
    Free tier users cannot access chat
    """
    db = SessionLocal()
    try:
        service = SubscriptionService(db)
        subscription = service.get_active_subscription(current_user.uid)

        if not subscription:
            # Free tier - check if free plan has chat
            plan = db.query(SubscriptionPlan).filter(
                SubscriptionPlan.name == 'Free'
            ).first()
        else:
            plan = service.get_plan(subscription.plan_id)

        if not plan or not plan.chat_support:
            raise HTTPException(
                status_code=status.HTTP_402_PAYMENT_REQUIRED,
                detail="Chat support is not available on your current plan. Please upgrade."
            )

        return {"plan": plan}
    finally:
        db.close()
