"""
Admin API Routes
Admin-only endpoints for platform management
"""

from fastapi import APIRouter, Depends, HTTPException, Query
from middleware.admin_guard import require_admin
from database import SessionLocal, Job, User
from database_extensions import Subscription, SubscriptionPlan, UsageTracking, ChatConversation, AdminUser, UserNote, UserTag, ChatMessage, WebhookEvent, AuditLog
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
            query = query.filter(Job.status == status)

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
                    if job.status in ["running", "pending"]:
                        job.status = "cancelled"
                        job.error_message = "Cancelled by admin"
                        results["success"].append(job_id)
                    else:
                        results["failed"].append({"job_id": job_id, "error": f"Cannot cancel job with status: {job.status}"})
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
                "status": job.status,
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
