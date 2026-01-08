"""
Admin API Routes
Admin-only endpoints for platform management
"""

from fastapi import APIRouter, Depends, HTTPException, Query
from middleware.admin_guard import require_admin
from database import SessionLocal, Job, User
from database_extensions import Subscription, SubscriptionPlan, UsageTracking, ChatConversation, AdminUser
from services.metrics_service import MetricsService
from services.audit_service import AuditService
from sqlalchemy import func
from datetime import datetime, timedelta
import os
from pathlib import Path
from typing import Optional

router = APIRouter(prefix="/admin", tags=["admin"])

@router.get("/dashboard/stats")
async def get_dashboard_stats(admin=Depends(require_admin)):
    """
    Get admin dashboard statistics
    """
    db = SessionLocal()
    try:
        # Count statistics
        total_users = db.query(User).count()
        active_subs = db.query(Subscription).filter(
            Subscription.status == 'active'
        ).count()

        total_jobs = db.query(Job).count()
        running_jobs = db.query(Job).filter(
            Job.status == 'running'
        ).count()

        # Current month revenue (approximate from subscriptions)
        subscriptions = db.query(Subscription, SubscriptionPlan).join(
            SubscriptionPlan
        ).filter(Subscription.status == 'active').all()

        monthly_revenue = sum(plan.price_cents for _, plan in subscriptions)

        # System metrics
        metrics = MetricsService.get_system_metrics()
        nextflow_count = MetricsService.get_nextflow_processes()

        # Unread chat count
        unread_chats = db.query(ChatConversation).filter(
            ChatConversation.status == 'open'
        ).count()

        return {
            "users": {
                "total": total_users,
                "active_subscriptions": active_subs
            },
            "jobs": {
                "total": total_jobs,
                "running": running_jobs
            },
            "revenue": {
                "monthly_recurring_cents": monthly_revenue,
                "monthly_recurring_usd": monthly_revenue / 100
            },
            "system": metrics,
            "nextflow_processes": nextflow_count,
            "unread_chats": unread_chats
        }
    finally:
        db.close()


@router.get("/jobs")
async def get_all_jobs(
    status: str = None,
    user_id: str = None,
    limit: int = 50,
    offset: int = 0,
    admin=Depends(require_admin)
):
    """
    Get all jobs across all users
    """
    db = SessionLocal()
    try:
        query = db.query(Job)

        if status:
            query = query.filter(Job.status == status)
        if user_id:
            query = query.filter(Job.user_id == user_id)

        total = query.count()
        jobs = query.order_by(Job.created_at.desc()).limit(limit).offset(offset).all()

        return {
            "jobs": jobs,
            "total": total,
            "limit": limit,
            "offset": offset
        }
    finally:
        db.close()


@router.get("/jobs/{job_id}/logs")
async def get_job_logs(job_id: str, admin=Depends(require_admin)):
    """
    Get Nextflow logs for a job
    """
    db = SessionLocal()
    try:
        job = db.query(Job).filter(Job.job_id == job_id).first()
        if not job:
            raise HTTPException(status_code=404, detail="Job not found")

        # Look for .nextflow.log in results directory
        results_dir = Path("./results") / job_id
        log_file = results_dir / ".nextflow.log"

        if log_file.exists():
            with open(log_file, 'r') as f:
                logs = f.read()
            return {"logs": logs, "path": str(log_file)}
        else:
            return {"logs": None, "error": "Log file not found"}
    finally:
        db.close()


@router.get("/users")
async def get_all_users(admin=Depends(require_admin)):
    """
    Get all users with subscription info
    """
    db = SessionLocal()
    try:
        users = db.query(User).all()

        user_data = []
        for user in users:
            subscription = db.query(Subscription).filter(
                Subscription.user_id == user.firebase_uid
            ).first()

            plan = None
            if subscription:
                plan = db.query(SubscriptionPlan).filter(
                    SubscriptionPlan.id == subscription.plan_id
                ).first()

            current_month = int(datetime.now().strftime('%Y%m'))
            usage = db.query(UsageTracking).filter(
                UsageTracking.user_id == user.firebase_uid,
                UsageTracking.month == current_month
            ).first()

            user_data.append({
                "uid": user.firebase_uid,
                "email": user.email,
                "created_at": user.created_at,
                "subscription": {
                    "plan": plan.name if plan else "Free",
                    "status": subscription.status if subscription else "none",
                } if subscription or plan else None,
                "usage": {
                    "jobs_executed": usage.jobs_executed if usage else 0,
                    "jobs_limit": usage.jobs_limit if usage else 2
                } if usage else None
            })

        return {"users": user_data, "total": len(user_data)}
    finally:
        db.close()


@router.post("/users/{user_uid}/usage")
async def update_user_usage(
    user_uid: str,
    jobs_limit: Optional[int] = Query(None),
    reset_count: bool = Query(False),
    admin=Depends(require_admin)
):
    """
    Update user's usage limits or reset usage count
    """
    # Store admin firebase_uid before session closes
    admin_firebase_uid = admin.firebase_uid

    db = SessionLocal()
    try:
        current_month = int(datetime.now().strftime('%Y%m'))

        # Get or create usage tracking for current month
        usage = db.query(UsageTracking).filter(
            UsageTracking.user_id == user_uid,
            UsageTracking.month == current_month
        ).first()

        if not usage:
            # Create new usage tracking entry
            usage = UsageTracking(
                user_id=user_uid,
                month=current_month,
                jobs_executed=0,
                jobs_limit=jobs_limit if jobs_limit is not None else 2
            )
            db.add(usage)
        else:
            # Update existing usage
            if jobs_limit is not None:
                usage.jobs_limit = jobs_limit
            if reset_count:
                usage.jobs_executed = 0

        db.commit()
        db.refresh(usage)

        # Log the action
        AuditService(db).log_action(
            action="update_user_usage",
            user_id=admin_firebase_uid,
            resource_type="user",
            resource_id=user_uid,
            metadata={
                "jobs_limit": usage.jobs_limit,
                "jobs_executed": usage.jobs_executed,
                "reset_count": reset_count
            }
        )

        return {
            "success": True,
            "usage": {
                "user_id": user_uid,
                "jobs_executed": usage.jobs_executed,
                "jobs_limit": usage.jobs_limit,
                "month": usage.month
            }
        }
    finally:
        db.close()


@router.post("/users/{user_uid}/subscription")
async def update_user_subscription(
    user_uid: str,
    plan_name: str = Query(..., description="Plan name: Free, Basic, Pro, or Enterprise"),
    admin=Depends(require_admin)
):
    """
    Update user subscription plan (admin only)
    """
    print(f"🔧 Admin subscription update request:")
    print(f"   User UID: {user_uid}")
    print(f"   Plan name: {plan_name}")
    print(f"   Admin: {admin.firebase_uid}")

    admin_firebase_uid = admin.firebase_uid
    db = SessionLocal()
    try:
        # Verify user exists
        user = db.query(User).filter(User.firebase_uid == user_uid).first()
        if not user:
            print(f"❌ User not found: {user_uid}")
            raise HTTPException(status_code=404, detail="User not found")

        print(f"✓ User found: {user.email}")

        # Get plan by name (case insensitive)
        plan = db.query(SubscriptionPlan).filter(
            func.lower(SubscriptionPlan.name) == plan_name.lower()
        ).first()

        if not plan:
            # List available plans for debugging
            all_plans = db.query(SubscriptionPlan).all()
            available_plans = [p.name for p in all_plans]
            print(f"❌ Plan '{plan_name}' not found")
            print(f"   Available plans: {available_plans}")
            raise HTTPException(
                status_code=404,
                detail=f"Plan '{plan_name}' not found. Available: {', '.join(available_plans) if available_plans else 'None - run database migration'}"
            )

        print(f"✓ Plan found: {plan.name} (ID: {plan.id})")

        # Check if user has existing subscription
        subscription = db.query(Subscription).filter(
            Subscription.user_id == user_uid
        ).first()

        current_period_start = datetime.utcnow()
        current_period_end = current_period_start + timedelta(days=30)

        if subscription:
            # Update existing subscription
            old_plan_id = subscription.plan_id
            subscription.plan_id = plan.id
            subscription.status = 'active'
            subscription.current_period_start = current_period_start
            subscription.current_period_end = current_period_end
            subscription.cancel_at_period_end = False
            subscription.updated_at = datetime.utcnow()
        else:
            # Create new subscription
            old_plan_id = None
            subscription = Subscription(
                user_id=user_uid,
                plan_id=plan.id,
                status='active',
                current_period_start=current_period_start,
                current_period_end=current_period_end,
                cancel_at_period_end=False
            )
            db.add(subscription)

        db.commit()
        db.refresh(subscription)

        # Update usage tracking with new plan limits
        current_month = int(datetime.now().strftime('%Y%m'))
        usage = db.query(UsageTracking).filter(
            UsageTracking.user_id == user_uid,
            UsageTracking.month == current_month
        ).first()

        if usage:
            usage.jobs_limit = plan.monthly_jobs_limit
        else:
            usage = UsageTracking(
                user_id=user_uid,
                month=current_month,
                jobs_executed=0,
                jobs_limit=plan.monthly_jobs_limit
            )
            db.add(usage)

        db.commit()

        # Log the action
        AuditService(db).log_action(
            action="admin_update_user_subscription",
            user_id=admin_firebase_uid,
            resource_type="subscription",
            resource_id=str(subscription.id),
            metadata={
                "user_uid": user_uid,
                "user_email": user.email,
                "old_plan_id": old_plan_id,
                "new_plan_id": plan.id,
                "new_plan_name": plan.name,
                "new_jobs_limit": plan.monthly_jobs_limit
            }
        )

        print(f"✅ Subscription updated successfully for {user.email}")

        return {
            "success": True,
            "subscription": {
                "id": subscription.id,
                "user_id": user_uid,
                "plan_name": plan.name,
                "plan_id": plan.id,
                "status": subscription.status,
                "monthly_jobs_limit": plan.monthly_jobs_limit,
                "price_usd": plan.price_cents / 100,
                "current_period_start": subscription.current_period_start,
                "current_period_end": subscription.current_period_end
            }
        }
    except HTTPException:
        # Re-raise HTTP exceptions (already handled)
        raise
    except Exception as e:
        # Catch any other database or unexpected errors
        print(f"❌ Unexpected error updating subscription: {type(e).__name__}: {str(e)}")
        import traceback
        traceback.print_exc()
        raise HTTPException(
            status_code=500,
            detail=f"Failed to update subscription: {str(e)}"
        )
    finally:
        db.close()


@router.get("/system/metrics")
async def get_system_metrics(admin=Depends(require_admin)):
    """
    Get real-time system metrics
    """
    metrics = MetricsService.get_system_metrics()
    nextflow_count = MetricsService.get_nextflow_processes()
    storage = MetricsService.get_storage_usage("./")

    return {
        "system": metrics,
        "nextflow_processes": nextflow_count,
        "storage": storage,
        "timestamp": datetime.utcnow().isoformat()
    }
