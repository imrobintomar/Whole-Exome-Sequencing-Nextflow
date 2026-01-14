"""
Admin API Routes
Admin-only endpoints for platform management
"""

from fastapi import APIRouter, Depends, HTTPException, Query
from middleware.admin_guard import require_admin
from database import SessionLocal, Job, User, JobStatus
from database_extensions import Subscription, SubscriptionPlan, UsageTracking, ChatConversation, AdminUser, UserNote, UserTag, ChatMessage, WebhookEvent, AuditLog
from services.metrics_service import MetricsService
from services.audit_service import AuditService
from services.email_service import get_email_service
from services.health_monitor import HealthMonitor
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
    search: str = None,
    date_from: str = None,
    date_to: str = None,
    limit: int = 50,
    offset: int = 0,
    admin=Depends(require_admin)
):
    """
    Get all jobs across all users with filtering
    """
    db = SessionLocal()
    try:
        query = db.query(Job)

        # Filter by status
        if status:
            # Convert string to enum
            try:
                status_enum = JobStatus(status.lower())
                query = query.filter(Job.status == status_enum)
            except ValueError:
                # Invalid status, skip filter
                pass

        # Filter by user
        if user_id:
            query = query.filter(Job.user_id == user_id)

        # Search by job_id or sample_name
        if search:
            query = query.filter(
                (Job.job_id.ilike(f"%{search}%")) |
                (Job.sample_name.ilike(f"%{search}%"))
            )

        # Filter by date range
        if date_from:
            query = query.filter(Job.created_at >= date_from)
        if date_to:
            query = query.filter(Job.created_at <= date_to)

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
                } if usage else None,
                "is_active": user.is_active if hasattr(user, 'is_active') else True,
                "is_banned": user.is_banned if hasattr(user, 'is_banned') else False,
                "ban_reason": user.ban_reason if hasattr(user, 'ban_reason') else None,
                "banned_at": user.banned_at if hasattr(user, 'banned_at') else None,
                "banned_by": user.banned_by if hasattr(user, 'banned_by') else None
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


@router.post("/jobs/bulk-action")
async def bulk_job_action(
    action: str,
    job_ids: list[str],
    admin=Depends(require_admin)
):
    """
    Perform bulk actions on multiple jobs
    Supported actions: cancel, delete
    """
    db = SessionLocal()
    try:
        results = {
            "success": [],
            "failed": []
        }

        for job_id in job_ids:
            try:
                job = db.query(Job).filter(Job.job_id == job_id).first()
                if not job:
                    results["failed"].append({"job_id": job_id, "error": "Not found"})
                    continue

                if action == "delete":
                    db.delete(job)
                    results["success"].append(job_id)

                elif action == "cancel":
                    if job.status in [JobStatus.RUNNING, JobStatus.PENDING]:
                        job.status = JobStatus.FAILED  # Mark as failed since we don't have CANCELLED status
                        job.error_message = "Cancelled by admin"
                        results["success"].append(job_id)
                    else:
                        results["failed"].append({"job_id": job_id, "error": f"Cannot cancel job with status: {job.status.value}"})
                else:
                    results["failed"].append({"job_id": job_id, "error": f"Unknown action: {action}"})

            except Exception as e:
                results["failed"].append({"job_id": job_id, "error": str(e)})

        db.commit()

        # Log the bulk action
        AuditService.log_action(
            admin.firebase_uid,
            f"bulk_{action}",
            {"job_ids": job_ids, "success": len(results["success"]), "failed": len(results["failed"])},
            f"Bulk {action} on {len(job_ids)} jobs"
        )

        return results
    finally:
        db.close()


@router.get("/users/export")
async def export_users(admin=Depends(require_admin)):
    """
    Export all users as CSV
    """
    from fastapi.responses import Response

    db = SessionLocal()
    try:
        users = db.query(User).all()

        # Create CSV header
        csv_data = "UID,Email,Username,Created At,Jobs Count,Active Subscription\n"

        for user in users:
            # Count jobs for each user
            jobs_count = db.query(Job).filter(Job.user_id == user.firebase_uid).count()

            # Check subscription
            subscription = db.query(Subscription).filter(
                Subscription.user_id == user.firebase_uid,
                Subscription.status == 'active'
            ).first()

            plan_name = "Free"
            if subscription:
                plan = db.query(SubscriptionPlan).filter(
                    SubscriptionPlan.id == subscription.plan_id
                ).first()
                if plan:
                    plan_name = plan.name

            # Escape commas in fields
            email = user.email.replace(',', ';') if user.email else ''
            username = user.username.replace(',', ';') if user.username else ''
            created_at = user.created_at.strftime('%Y-%m-%d %H:%M:%S') if user.created_at else ''

            csv_data += f"{user.firebase_uid},{email},{username},{created_at},{jobs_count},{plan_name}\n"

        return Response(
            content=csv_data,
            media_type="text/csv",
            headers={
                "Content-Disposition": f"attachment; filename=users_export_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
            }
        )
    finally:
        db.close()


@router.post("/users/{user_uid}/ban")
async def ban_user(
    user_uid: str,
    reason: str = Query(..., description="Reason for banning"),
    admin=Depends(require_admin)
):
    """Ban a user account"""
    db = SessionLocal()
    try:
        user = db.query(User).filter(User.firebase_uid == user_uid).first()
        if not user:
            raise HTTPException(status_code=404, detail="User not found")

        user.is_banned = True
        user.ban_reason = reason
        user.banned_at = datetime.utcnow()
        user.banned_by = admin.firebase_uid
        db.commit()

        AuditService.log_action(admin.firebase_uid, "ban_user", {"user_uid": user_uid, "reason": reason}, f"Banned user {user.email}")

        return {"success": True, "message": f"User {user.email} has been banned", "user_uid": user_uid}
    finally:
        db.close()


@router.post("/users/{user_uid}/unban")
async def unban_user(user_uid: str, admin=Depends(require_admin)):
    """Unban a user account"""
    db = SessionLocal()
    try:
        user = db.query(User).filter(User.firebase_uid == user_uid).first()
        if not user:
            raise HTTPException(status_code=404, detail="User not found")

        user.is_banned = False
        user.ban_reason = None
        user.banned_at = None
        user.banned_by = None
        db.commit()

        AuditService.log_action(admin.firebase_uid, "unban_user", {"user_uid": user_uid}, f"Unbanned user {user.email}")

        return {"success": True, "message": f"User {user.email} has been unbanned"}
    finally:
        db.close()


@router.post("/users/{user_uid}/suspend")
async def suspend_user(user_uid: str, admin=Depends(require_admin)):
    """Suspend a user account (set is_active to False)"""
    db = SessionLocal()
    try:
        user = db.query(User).filter(User.firebase_uid == user_uid).first()
        if not user:
            raise HTTPException(status_code=404, detail="User not found")

        user.is_active = False
        db.commit()

        AuditService.log_action(admin.firebase_uid, "suspend_user", {"user_uid": user_uid}, f"Suspended user {user.email}")

        return {"success": True, "message": f"User {user.email} has been suspended"}
    finally:
        db.close()


@router.post("/users/{user_uid}/activate")
async def activate_user(user_uid: str, admin=Depends(require_admin)):
    """Activate a suspended user account (set is_active to True)"""
    db = SessionLocal()
    try:
        user = db.query(User).filter(User.firebase_uid == user_uid).first()
        if not user:
            raise HTTPException(status_code=404, detail="User not found")

        user.is_active = True
        db.commit()

        AuditService.log_action(admin.firebase_uid, "activate_user", {"user_uid": user_uid}, f"Activated user {user.email}")

        return {"success": True, "message": f"User {user.email} has been activated"}
    finally:
        db.close()


@router.get("/users/{user_uid}/details")
async def get_user_details(user_uid: str, admin=Depends(require_admin)):
    """
    Get comprehensive user details including:
    - User profile and subscription
    - Full job history
    - Payment history (Stripe events)
    - Support ticket history
    - Activity timeline (audit logs)
    - Notes and tags
    """
    db = SessionLocal()
    try:
        # Get user
        user = db.query(User).filter(User.firebase_uid == user_uid).first()
        if not user:
            raise HTTPException(status_code=404, detail="User not found")

        # Get subscription details
        subscription = db.query(Subscription).filter(
            Subscription.user_id == user_uid
        ).first()

        plan = None
        if subscription:
            plan = db.query(SubscriptionPlan).filter(
                SubscriptionPlan.id == subscription.plan_id
            ).first()

        # Get usage tracking
        current_month = int(datetime.now().strftime('%Y%m'))
        usage = db.query(UsageTracking).filter(
            UsageTracking.user_id == user_uid,
            UsageTracking.month == current_month
        ).first()

        # Get all jobs
        jobs = db.query(Job).filter(Job.user_id == user_uid).order_by(Job.created_at.desc()).all()

        # Get payment history (webhook events for this customer)
        payment_history = []
        if subscription and subscription.stripe_customer_id:
            webhooks = db.query(WebhookEvent).filter(
                WebhookEvent.event_type.in_([
                    'payment_intent.succeeded',
                    'payment_intent.payment_failed',
                    'invoice.paid',
                    'invoice.payment_failed',
                    'charge.succeeded',
                    'charge.failed'
                ])
            ).order_by(WebhookEvent.created_at.desc()).limit(50).all()

            # Filter webhooks for this customer (by checking payload)
            for webhook in webhooks:
                import json
                try:
                    payload = json.loads(webhook.payload_json)
                    customer_id = None

                    # Extract customer_id from different event types
                    if 'data' in payload and 'object' in payload['data']:
                        obj = payload['data']['object']
                        customer_id = obj.get('customer')

                    if customer_id == subscription.stripe_customer_id:
                        payment_history.append({
                            "id": webhook.id,
                            "event_type": webhook.event_type,
                            "created_at": webhook.created_at,
                            "processed": webhook.processed,
                            "amount": obj.get('amount', 0) / 100 if 'amount' in obj else None,
                            "currency": obj.get('currency', 'usd')
                        })
                except:
                    pass

        # Get support tickets
        conversations = db.query(ChatConversation).filter(
            ChatConversation.user_id == user_uid
        ).order_by(ChatConversation.created_at.desc()).all()

        support_tickets = []
        for conv in conversations:
            message_count = db.query(ChatMessage).filter(
                ChatMessage.conversation_id == conv.id
            ).count()

            support_tickets.append({
                "id": conv.id,
                "subject": conv.subject,
                "status": conv.status,
                "job_id": conv.job_id,
                "message_count": message_count,
                "last_message_at": conv.last_message_at,
                "created_at": conv.created_at
            })

        # Get activity timeline (audit logs)
        audit_logs = db.query(AuditLog).filter(
            (AuditLog.user_id == user_uid) |
            (AuditLog.resource_id == user_uid)
        ).order_by(AuditLog.created_at.desc()).limit(100).all()

        activity_timeline = [{
            "id": log.id,
            "action": log.action,
            "resource_type": log.resource_type,
            "resource_id": log.resource_id,
            "created_at": log.created_at,
            "metadata": log.metadata_json
        } for log in audit_logs]

        # Get notes
        notes = db.query(UserNote).filter(
            UserNote.user_id == user_uid
        ).order_by(UserNote.created_at.desc()).all()

        user_notes = [{
            "id": note.id,
            "note_text": note.note_text,
            "admin_id": note.admin_id,
            "created_at": note.created_at,
            "updated_at": note.updated_at
        } for note in notes]

        # Get tags
        tags = db.query(UserTag).filter(
            UserTag.user_id == user_uid
        ).all()

        user_tags = [{
            "id": tag.id,
            "tag_name": tag.tag_name,
            "color": tag.color,
            "created_at": tag.created_at,
            "created_by": tag.created_by
        } for tag in tags]

        # Compile response
        return {
            "user": {
                "uid": user.firebase_uid,
                "email": user.email,
                "username": user.username,
                "created_at": user.created_at,
                "is_active": user.is_active if hasattr(user, 'is_active') else True,
                "is_banned": user.is_banned if hasattr(user, 'is_banned') else False,
                "ban_reason": user.ban_reason if hasattr(user, 'ban_reason') else None,
                "banned_at": user.banned_at if hasattr(user, 'banned_at') else None,
            },
            "subscription": {
                "plan": plan.name if plan else "Free",
                "status": subscription.status if subscription else "none",
                "stripe_customer_id": subscription.stripe_customer_id if subscription else None,
                "current_period_start": subscription.current_period_start if subscription else None,
                "current_period_end": subscription.current_period_end if subscription else None,
                "price_cents": plan.price_cents if plan else 0,
                "monthly_jobs_limit": plan.monthly_jobs_limit if plan else 2,
            },
            "usage": {
                "jobs_executed": usage.jobs_executed if usage else 0,
                "jobs_limit": usage.jobs_limit if usage else 2,
                "month": usage.month if usage else current_month
            },
            "jobs": [{
                "job_id": job.job_id,
                "sample_name": job.sample_name,
                "status": job.status.value if hasattr(job.status, 'value') else str(job.status),
                "current_step": job.current_step,
                "error_message": job.error_message,
                "created_at": job.created_at,
                "updated_at": job.updated_at,
                "completed_at": job.completed_at if hasattr(job, 'completed_at') else None
            } for job in jobs],
            "payment_history": payment_history,
            "support_tickets": support_tickets,
            "activity_timeline": activity_timeline,
            "notes": user_notes,
            "tags": user_tags
        }
    finally:
        db.close()


@router.post("/users/{user_uid}/notes")
async def add_user_note(
    user_uid: str,
    note_text: str = Query(..., description="Note content"),
    admin=Depends(require_admin)
):
    """Add a note to a user's profile"""
    db = SessionLocal()
    try:
        # Verify user exists
        user = db.query(User).filter(User.firebase_uid == user_uid).first()
        if not user:
            raise HTTPException(status_code=404, detail="User not found")

        # Create note
        note = UserNote(
            user_id=user_uid,
            admin_id=admin.firebase_uid,
            note_text=note_text
        )
        db.add(note)
        db.commit()
        db.refresh(note)

        AuditService(db).log_action(
            action="add_user_note",
            user_id=admin.firebase_uid,
            resource_type="user",
            resource_id=user_uid,
            metadata={"note_preview": note_text[:50]}
        )

        return {
            "success": True,
            "note": {
                "id": note.id,
                "note_text": note.note_text,
                "admin_id": note.admin_id,
                "created_at": note.created_at
            }
        }
    finally:
        db.close()


@router.post("/users/{user_uid}/tags")
async def add_user_tag(
    user_uid: str,
    tag_name: str = Query(..., description="Tag name"),
    color: str = Query("blue", description="Tag color"),
    admin=Depends(require_admin)
):
    """Add a tag to a user's profile"""
    db = SessionLocal()
    try:
        # Verify user exists
        user = db.query(User).filter(User.firebase_uid == user_uid).first()
        if not user:
            raise HTTPException(status_code=404, detail="User not found")

        # Check if tag already exists
        existing_tag = db.query(UserTag).filter(
            UserTag.user_id == user_uid,
            UserTag.tag_name == tag_name
        ).first()

        if existing_tag:
            raise HTTPException(status_code=400, detail="Tag already exists for this user")

        # Create tag
        tag = UserTag(
            user_id=user_uid,
            tag_name=tag_name,
            color=color,
            created_by=admin.firebase_uid
        )
        db.add(tag)
        db.commit()
        db.refresh(tag)

        return {
            "success": True,
            "tag": {
                "id": tag.id,
                "tag_name": tag.tag_name,
                "color": tag.color,
                "created_at": tag.created_at
            }
        }
    finally:
        db.close()


@router.delete("/users/{user_uid}/tags/{tag_id}")
async def remove_user_tag(
    user_uid: str,
    tag_id: int,
    admin=Depends(require_admin)
):
    """Remove a tag from a user's profile"""
    db = SessionLocal()
    try:
        tag = db.query(UserTag).filter(
            UserTag.id == tag_id,
            UserTag.user_id == user_uid
        ).first()

        if not tag:
            raise HTTPException(status_code=404, detail="Tag not found")

        db.delete(tag)
        db.commit()

        return {"success": True, "message": "Tag removed"}
    finally:
        db.close()


# ============================================================================
# EMAIL NOTIFICATION ENDPOINTS
# ============================================================================

@router.post("/email/test")
async def test_email_connection(admin=Depends(require_admin)):
    """
    Test SMTP email connection
    """
    email_service = get_email_service()
    result = email_service.test_connection()

    AuditService.log_action(
        admin.firebase_uid,
        "test_email",
        result,
        "Tested email connection"
    )

    return result


@router.post("/email/send-custom")
async def send_custom_email(
    user_email: str = Query(..., description="Recipient email"),
    subject: str = Query(..., description="Email subject"),
    message: str = Query(..., description="Email message"),
    user_name: str = Query(None, description="Optional user name"),
    admin=Depends(require_admin)
):
    """
    Send custom notification email to a user
    """
    try:
        email_service = get_email_service()
        success = email_service.send_custom_notification(
            user_email=user_email,
            subject=subject,
            message=message,
            user_name=user_name
        )

        # Log the action
        try:
            AuditService.log_action(
                admin.firebase_uid,
                "send_custom_email",
                {"user_email": user_email, "subject": subject, "success": success},
                f"Sent custom email to {user_email}"
            )
        except Exception as audit_error:
            print(f"⚠️ Failed to log audit: {audit_error}")

        if success:
            return {"success": True, "message": f"Email sent to {user_email}"}
        else:
            raise HTTPException(status_code=500, detail="Failed to send email - check SMTP settings and recipient address")

    except HTTPException:
        raise
    except Exception as e:
        print(f"❌ Error in send_custom_email: {str(e)}")
        import traceback
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=f"Failed to send email: {str(e)}")


@router.post("/email/payment-reminder")
async def send_payment_reminder(
    user_uid: str = Query(..., description="User UID"),
    amount: float = Query(..., description="Payment amount"),
    due_date: str = Query(..., description="Payment due date"),
    invoice_url: str = Query(None, description="Optional invoice URL"),
    admin=Depends(require_admin)
):
    """
    Send payment due reminder to a user
    """
    db = SessionLocal()
    try:
        # Get user details
        user = db.query(User).filter(User.firebase_uid == user_uid).first()
        if not user:
            raise HTTPException(status_code=404, detail="User not found")

        email_service = get_email_service()
        success = email_service.send_payment_due_alert(
            user_email=user.email,
            user_name=user.username or user.email.split('@')[0],
            amount=amount,
            due_date=due_date,
            invoice_url=invoice_url
        )

        AuditService.log_action(
            admin.firebase_uid,
            "send_payment_reminder",
            {"user_uid": user_uid, "amount": amount, "success": success},
            f"Sent payment reminder to {user.email}"
        )

        if success:
            return {"success": True, "message": f"Payment reminder sent to {user.email}"}
        else:
            raise HTTPException(status_code=500, detail="Failed to send email")
    finally:
        db.close()


@router.post("/email/subscription-expiry")
async def send_subscription_expiry(
    user_uid: str = Query(..., description="User UID"),
    admin=Depends(require_admin)
):
    """
    Send subscription expiry notification to a user
    """
    db = SessionLocal()
    try:
        # Get user details
        user = db.query(User).filter(User.firebase_uid == user_uid).first()
        if not user:
            raise HTTPException(status_code=404, detail="User not found")

        # Get subscription details
        subscription = db.query(Subscription).filter(
            Subscription.user_id == user_uid
        ).first()

        if not subscription:
            raise HTTPException(status_code=404, detail="No active subscription found")

        plan = db.query(SubscriptionPlan).filter(
            SubscriptionPlan.id == subscription.plan_id
        ).first()

        email_service = get_email_service()
        success = email_service.send_subscription_expiry_alert(
            user_email=user.email,
            user_name=user.username or user.email.split('@')[0],
            expiry_date=subscription.current_period_end.strftime('%Y-%m-%d') if subscription.current_period_end else "Unknown",
            plan_name=plan.name if plan else "Unknown"
        )

        AuditService.log_action(
            admin.firebase_uid,
            "send_subscription_expiry",
            {"user_uid": user_uid, "success": success},
            f"Sent subscription expiry notice to {user.email}"
        )

        if success:
            return {"success": True, "message": f"Subscription expiry notice sent to {user.email}"}
        else:
            raise HTTPException(status_code=500, detail="Failed to send email")
    finally:
        db.close()


@router.get("/email/health-alerts")
async def get_health_alerts(
    limit: int = Query(50, description="Max number of alerts"),
    severity: str = Query(None, description="Filter by severity (warning/critical)"),
    admin=Depends(require_admin)
):
    """
    Get health alert history
    """
    alerts = HealthMonitor.get_alert_history(limit=limit, severity=severity)
    return {
        "alerts": alerts,
        "count": len(alerts)
    }


@router.post("/email/health-alerts/{alert_id}/resolve")
async def resolve_health_alert(
    alert_id: int,
    admin=Depends(require_admin)
):
    """
    Mark a health alert as resolved
    """
    success = HealthMonitor.resolve_alert(alert_id)

    if success:
        AuditService.log_action(
            admin.firebase_uid,
            "resolve_health_alert",
            {"alert_id": alert_id},
            f"Resolved health alert #{alert_id}"
        )
        return {"success": True, "message": f"Alert #{alert_id} resolved"}
    else:
        raise HTTPException(status_code=404, detail="Alert not found")


@router.post("/email/health-alerts/test")
async def test_health_alert_email(admin=Depends(require_admin)):
    """
    Send a test health alert email
    """
    email_service = get_email_service()
    success = email_service.send_health_alert(
        alert_type="test",
        severity="warning",
        message="This is a test health alert email",
        value=85.5,
        threshold=80.0,
        metrics={"cpu": 45.2, "memory": 62.8, "disk": 55.1}
    )

    AuditService.log_action(
        admin.firebase_uid,
        "test_health_alert_email",
        {"success": success},
        "Sent test health alert email"
    )

    if success:
        return {"success": True, "message": "Test health alert email sent"}
    else:
        raise HTTPException(status_code=500, detail="Failed to send test email")


@router.get("/email/settings")
async def get_email_settings(admin=Depends(require_admin)):
    """
    Get current email notification settings
    """
    return {
        "smtp_host": os.getenv("SMTP_HOST", "smtp.hostinger.com"),
        "smtp_port": int(os.getenv("SMTP_PORT", "465")),
        "smtp_user": os.getenv("SMTP_USER", ""),
        "admin_email": os.getenv("ADMIN_EMAIL", ""),
        "health_alerts_enabled": os.getenv("HEALTH_EMAIL_ALERTS_ENABLED", "true").lower() == "true",
        "health_thresholds": {
            "cpu": int(os.getenv("HEALTH_CPU_THRESHOLD", "90")),
            "memory": int(os.getenv("HEALTH_MEMORY_THRESHOLD", "90")),
            "disk": int(os.getenv("HEALTH_DISK_THRESHOLD", "90"))
        }
    }


# ============================================================================
# ANALYTICS ENDPOINTS
# ============================================================================

@router.get("/analytics/summary")
async def get_analytics_summary(
    period: str = Query("month", description="Period: day, week, month, year"),
    admin=Depends(require_admin)
):
    """
    Get analytics summary for dashboard
    """
    db = SessionLocal()
    try:
        from datetime import datetime, timedelta

        # Calculate date range based on period
        now = datetime.utcnow()
        if period == "day":
            start_date = now - timedelta(days=1)
        elif period == "week":
            start_date = now - timedelta(weeks=1)
        elif period == "year":
            start_date = now - timedelta(days=365)
        else:  # month
            start_date = now - timedelta(days=30)

        # User analytics
        total_users = db.query(User).count()
        new_users_period = db.query(User).filter(User.created_at >= start_date).count()

        # Subscription analytics
        active_subs = db.query(Subscription).filter(Subscription.status == 'active').count()

        subscriptions = db.query(Subscription, SubscriptionPlan).join(
            SubscriptionPlan
        ).filter(Subscription.status == 'active').all()

        mrr = sum(plan.price_cents for _, plan in subscriptions) / 100

        # Job analytics
        total_jobs = db.query(Job).count()
        jobs_period = db.query(Job).filter(Job.created_at >= start_date).count()

        completed_jobs = db.query(Job).filter(
            Job.status == 'completed',
            Job.created_at >= start_date
        ).count()

        failed_jobs = db.query(Job).filter(
            Job.status == 'failed',
            Job.created_at >= start_date
        ).count()

        success_rate = (completed_jobs / jobs_period * 100) if jobs_period > 0 else 0

        return {
            "users": {
                "total": total_users,
                "new_this_period": new_users_period,
                "active_subscriptions": active_subs
            },
            "revenue": {
                "mrr": mrr,
                "total_subscriptions": active_subs
            },
            "jobs": {
                "total": total_jobs,
                "this_period": jobs_period,
                "completed": completed_jobs,
                "failed": failed_jobs,
                "success_rate": success_rate
            },
            "period": period,
            "start_date": start_date.isoformat(),
            "end_date": now.isoformat()
        }
    finally:
        db.close()


@router.get("/analytics/revenue")
async def get_revenue_analytics(
    from_date: str = Query(None, description="Start date (YYYY-MM-DD)"),
    to_date: str = Query(None, description="End date (YYYY-MM-DD)"),
    admin=Depends(require_admin)
):
    """
    Get revenue analytics over time
    """
    db = SessionLocal()
    try:
        from datetime import datetime, timedelta
        from collections import defaultdict

        # Parse dates or use defaults (last 12 months)
        if from_date and to_date:
            start = datetime.fromisoformat(from_date)
            end = datetime.fromisoformat(to_date)
        else:
            end = datetime.utcnow()
            start = end - timedelta(days=365)

        # Get all subscriptions
        subscriptions = db.query(Subscription, SubscriptionPlan).join(
            SubscriptionPlan
        ).filter(
            Subscription.created_at >= start,
            Subscription.created_at <= end
        ).all()

        # Group by month
        monthly_revenue = defaultdict(float)
        revenue_by_plan = defaultdict(float)

        for sub, plan in subscriptions:
            month_key = sub.created_at.strftime('%Y-%m')
            monthly_revenue[month_key] += plan.price_cents / 100
            revenue_by_plan[plan.name] += plan.price_cents / 100

        # Format for frontend
        revenue_over_time = [
            {"month": month, "revenue": revenue}
            for month, revenue in sorted(monthly_revenue.items())
        ]

        revenue_by_plan_data = [
            {"plan": plan, "revenue": revenue}
            for plan, revenue in revenue_by_plan.items()
        ]

        # Current MRR
        active_subs = db.query(Subscription, SubscriptionPlan).join(
            SubscriptionPlan
        ).filter(Subscription.status == 'active').all()

        current_mrr = sum(plan.price_cents for _, plan in active_subs) / 100

        # Calculate ARPU (Average Revenue Per User)
        total_users = db.query(User).count()
        arpu = current_mrr / total_users if total_users > 0 else 0

        return {
            "revenue_over_time": revenue_over_time,
            "revenue_by_plan": revenue_by_plan_data,
            "current_mrr": current_mrr,
            "arpu": arpu,
            "total_revenue": sum(monthly_revenue.values()),
            "from_date": start.isoformat(),
            "to_date": end.isoformat()
        }
    finally:
        db.close()


@router.get("/analytics/users")
async def get_user_analytics(
    from_date: str = Query(None, description="Start date (YYYY-MM-DD)"),
    to_date: str = Query(None, description="End date (YYYY-MM-DD)"),
    admin=Depends(require_admin)
):
    """
    Get user growth analytics
    """
    db = SessionLocal()
    try:
        from datetime import datetime, timedelta
        from collections import defaultdict

        # Parse dates or use defaults (last 12 months)
        if from_date and to_date:
            start = datetime.fromisoformat(from_date)
            end = datetime.fromisoformat(to_date)
        else:
            end = datetime.utcnow()
            start = end - timedelta(days=365)

        # Get all users
        users = db.query(User).filter(
            User.created_at >= start,
            User.created_at <= end
        ).all()

        # Group by month
        monthly_users = defaultdict(int)
        for user in users:
            month_key = user.created_at.strftime('%Y-%m')
            monthly_users[month_key] += 1

        # Cumulative growth
        user_growth = []
        cumulative = 0
        for month in sorted(monthly_users.keys()):
            cumulative += monthly_users[month]
            user_growth.append({
                "month": month,
                "new_users": monthly_users[month],
                "total_users": cumulative
            })

        # Users by plan
        users_by_plan = db.query(
            SubscriptionPlan.name,
            func.count(Subscription.id).label('count')
        ).join(
            Subscription, SubscriptionPlan.id == Subscription.plan_id
        ).filter(
            Subscription.status == 'active'
        ).group_by(SubscriptionPlan.name).all()

        free_users = db.query(User).count() - db.query(Subscription).filter(
            Subscription.status == 'active'
        ).count()

        users_by_plan_data = [{"plan": name, "count": count} for name, count in users_by_plan]
        users_by_plan_data.append({"plan": "Free", "count": free_users})

        return {
            "user_growth": user_growth,
            "users_by_plan": users_by_plan_data,
            "total_users": db.query(User).count(),
            "from_date": start.isoformat(),
            "to_date": end.isoformat()
        }
    finally:
        db.close()


@router.get("/analytics/jobs")
async def get_job_analytics(
    from_date: str = Query(None, description="Start date (YYYY-MM-DD)"),
    to_date: str = Query(None, description="End date (YYYY-MM-DD)"),
    admin=Depends(require_admin)
):
    """
    Get job analytics
    """
    db = SessionLocal()
    try:
        from datetime import datetime, timedelta
        from collections import defaultdict

        # Parse dates or use defaults (last 30 days)
        if from_date and to_date:
            start = datetime.fromisoformat(from_date)
            end = datetime.fromisoformat(to_date)
        else:
            end = datetime.utcnow()
            start = end - timedelta(days=30)

        # Get all jobs in period
        jobs = db.query(Job).filter(
            Job.created_at >= start,
            Job.created_at <= end
        ).all()

        # Group by day and status
        jobs_by_day = defaultdict(lambda: {"completed": 0, "failed": 0, "running": 0, "pending": 0})
        jobs_by_status = defaultdict(int)

        for job in jobs:
            day_key = job.created_at.strftime('%Y-%m-%d')
            status_str = job.status.value if hasattr(job.status, 'value') else str(job.status)
            jobs_by_day[day_key][status_str] += 1
            jobs_by_status[status_str] += 1

        # Format for frontend
        jobs_over_time = [
            {
                "date": day,
                "completed": data["completed"],
                "failed": data["failed"],
                "running": data["running"],
                "pending": data["pending"]
            }
            for day, data in sorted(jobs_by_day.items())
        ]

        jobs_by_status_data = [
            {"status": status, "count": count}
            for status, count in jobs_by_status.items()
        ]

        # Calculate success rate
        total_finished = jobs_by_status["completed"] + jobs_by_status["failed"]
        success_rate = (jobs_by_status["completed"] / total_finished * 100) if total_finished > 0 else 0

        # Calculate average job duration (if available)
        completed_jobs = [j for j in jobs if (j.status.value if hasattr(j.status, 'value') else str(j.status)) == 'completed' and hasattr(j, 'completed_at') and j.completed_at]
        avg_duration = 0
        if completed_jobs:
            durations = [(j.completed_at - j.created_at).total_seconds() / 3600 for j in completed_jobs]
            avg_duration = sum(durations) / len(durations)

        return {
            "jobs_over_time": jobs_over_time,
            "jobs_by_status": jobs_by_status_data,
            "success_rate": success_rate,
            "avg_duration_hours": avg_duration,
            "total_jobs": len(jobs),
            "from_date": start.isoformat(),
            "to_date": end.isoformat()
        }
    finally:
        db.close()
