"""
Billing API Routes
Subscription and payment management endpoints
"""

from fastapi import APIRouter, Depends, HTTPException, Request, Header
from firebase_auth import get_current_user
from database import SessionLocal, User
from database_extensions import SubscriptionPlan, Subscription, UsageTracking
from services.stripe_service import StripeService
from services.audit_service import AuditService
from datetime import datetime
import os
import json

router = APIRouter(prefix="/billing", tags=["billing"])


@router.get("/plans")
async def get_subscription_plans():
    """
    Get all available subscription plans (public endpoint)
    """
    db = SessionLocal()
    try:
        plans = db.query(SubscriptionPlan).filter(
            SubscriptionPlan.active == True
        ).all()

        return {
            "plans": [
                {
                    "id": plan.id,
                    "name": plan.name,
                    "price_cents": plan.price_cents,
                    "price_usd": plan.price_cents / 100,
                    "monthly_jobs_limit": plan.monthly_jobs_limit,
                    "chat_support": plan.chat_support,
                    "features": json.loads(plan.features_json) if plan.features_json else []
                }
                for plan in plans
            ]
        }
    finally:
        db.close()


@router.get("/subscription")
async def get_my_subscription(current_user=Depends(get_current_user)):
    """
    Get current user's subscription details
    """
    db = SessionLocal()
    try:
        subscription = db.query(Subscription).filter(
            Subscription.user_id == current_user.uid
        ).first()

        if not subscription:
            # Return free tier info
            free_plan = db.query(SubscriptionPlan).filter(
                SubscriptionPlan.name == "Free"
            ).first()

            current_month = int(datetime.now().strftime('%Y%m'))
            usage = db.query(UsageTracking).filter(
                UsageTracking.user_id == current_user.uid,
                UsageTracking.month == current_month
            ).first()

            return {
                "subscription": None,
                "plan": {
                    "name": "Free",
                    "monthly_jobs_limit": free_plan.monthly_jobs_limit if free_plan else 2,
                    "chat_support": False,
                    "price_cents": 0
                },
                "usage": {
                    "jobs_executed": usage.jobs_executed if usage else 0,
                    "jobs_limit": free_plan.monthly_jobs_limit if free_plan else 2,
                    "month": current_month
                }
            }

        # Get plan details
        plan = db.query(SubscriptionPlan).filter(
            SubscriptionPlan.id == subscription.plan_id
        ).first()

        # Get current month usage
        current_month = int(datetime.now().strftime('%Y%m'))
        usage = db.query(UsageTracking).filter(
            UsageTracking.user_id == current_user.uid,
            UsageTracking.month == current_month
        ).first()

        return {
            "subscription": {
                "status": subscription.status,
                "current_period_start": subscription.current_period_start,
                "current_period_end": subscription.current_period_end,
                "cancel_at_period_end": subscription.cancel_at_period_end,
                "stripe_customer_id": subscription.stripe_customer_id,
                "stripe_subscription_id": subscription.stripe_subscription_id
            },
            "plan": {
                "name": plan.name,
                "monthly_jobs_limit": plan.monthly_jobs_limit,
                "chat_support": plan.chat_support,
                "price_cents": plan.price_cents,
                "features": json.loads(plan.features_json) if plan.features_json else []
            },
            "usage": {
                "jobs_executed": usage.jobs_executed if usage else 0,
                "jobs_limit": usage.jobs_limit if usage else plan.monthly_jobs_limit,
                "month": current_month
            }
        }
    finally:
        db.close()


@router.post("/checkout")
async def create_checkout_session(
    plan_id: int,
    current_user=Depends(get_current_user)
):
    """
    Create Stripe checkout session for subscription
    """
    db = SessionLocal()
    try:
        # Verify plan exists
        plan = db.query(SubscriptionPlan).filter(
            SubscriptionPlan.id == plan_id,
            SubscriptionPlan.active == True
        ).first()

        if not plan:
            raise HTTPException(status_code=404, detail="Plan not found")

        if not plan.stripe_price_id:
            raise HTTPException(
                status_code=400,
                detail="Stripe price ID not configured for this plan"
            )

        # Check if user already has active subscription
        existing = db.query(Subscription).filter(
            Subscription.user_id == current_user.uid,
            Subscription.status.in_(['active', 'trialing'])
        ).first()

        if existing:
            raise HTTPException(
                status_code=400,
                detail="User already has an active subscription. Please cancel first."
            )

        # Get frontend URL from ENV
        frontend_url = os.getenv("FRONTEND_URL", "http://localhost:3000")
        success_url = f"{frontend_url}/dashboard?checkout=success"
        cancel_url = f"{frontend_url}/pricing?checkout=cancelled"

        # Create Stripe checkout session
        stripe_service = StripeService()
        checkout_url = stripe_service.create_checkout_session(
            price_id=plan.stripe_price_id,
            customer_email=current_user.email,
            success_url=success_url,
            cancel_url=cancel_url,
            metadata={
                "user_id": current_user.uid,
                "plan_id": str(plan_id),
                "plan_name": plan.name
            }
        )

        # Audit log
        audit = AuditService(db)
        audit.log_action(
            action="checkout_initiated",
            user_id=current_user.uid,
            resource_type="subscription",
            metadata={"plan_id": plan_id, "plan_name": plan.name}
        )

        return {
            "checkout_url": checkout_url,
            "plan": {
                "name": plan.name,
                "price_cents": plan.price_cents
            }
        }
    finally:
        db.close()


@router.post("/portal")
async def create_customer_portal_session(current_user=Depends(get_current_user)):
    """
    Create Stripe customer portal session for managing subscription
    """
    db = SessionLocal()
    try:
        subscription = db.query(Subscription).filter(
            Subscription.user_id == current_user.uid
        ).first()

        if not subscription or not subscription.stripe_customer_id:
            raise HTTPException(
                status_code=404,
                detail="No subscription found. Please subscribe first."
            )

        frontend_url = os.getenv("FRONTEND_URL", "http://localhost:3000")
        return_url = f"{frontend_url}/dashboard/billing"

        stripe_service = StripeService()
        portal_url = stripe_service.create_customer_portal_session(
            customer_id=subscription.stripe_customer_id,
            return_url=return_url
        )

        # Audit log
        audit = AuditService(db)
        audit.log_action(
            action="customer_portal_accessed",
            user_id=current_user.uid,
            resource_type="subscription"
        )

        return {"portal_url": portal_url}
    finally:
        db.close()


@router.get("/usage")
async def get_usage_stats(current_user=Depends(get_current_user)):
    """
    Get detailed usage statistics for current month
    """
    db = SessionLocal()
    try:
        current_month = int(datetime.now().strftime('%Y%m'))

        usage = db.query(UsageTracking).filter(
            UsageTracking.user_id == current_user.uid,
            UsageTracking.month == current_month
        ).first()

        subscription = db.query(Subscription).filter(
            Subscription.user_id == current_user.uid
        ).first()

        # Get plan
        if subscription:
            plan = db.query(SubscriptionPlan).filter(
                SubscriptionPlan.id == subscription.plan_id
            ).first()
        else:
            plan = db.query(SubscriptionPlan).filter(
                SubscriptionPlan.name == "Free"
            ).first()

        jobs_executed = usage.jobs_executed if usage else 0
        jobs_limit = plan.monthly_jobs_limit if plan else 2

        return {
            "month": current_month,
            "jobs_executed": jobs_executed,
            "jobs_limit": jobs_limit,
            "jobs_remaining": max(0, jobs_limit - jobs_executed),
            "usage_percent": (jobs_executed / jobs_limit * 100) if jobs_limit > 0 else 0,
            "plan_name": plan.name if plan else "Free",
            "chat_support_enabled": plan.chat_support if plan else False
        }
    finally:
        db.close()


@router.post("/webhooks/stripe")
async def stripe_webhook(
    request: Request,
    stripe_signature: str = Header(None, alias="stripe-signature")
):
    """
    Handle Stripe webhook events
    IMPORTANT: This endpoint must be publicly accessible (no auth)
    """
    db = SessionLocal()
    try:
        # Get raw body
        payload = await request.body()

        # Verify webhook signature
        stripe_service = StripeService()
        event = stripe_service.verify_webhook(
            payload=payload,
            signature=stripe_signature
        )

        if not event:
            raise HTTPException(status_code=400, detail="Invalid signature")

        # Check for duplicate webhooks
        from database_extensions import WebhookEvent
        existing = db.query(WebhookEvent).filter(
            WebhookEvent.stripe_event_id == event['id']
        ).first()

        if existing:
            return {"status": "duplicate", "event_id": event['id']}

        # Log webhook event
        webhook_log = WebhookEvent(
            stripe_event_id=event['id'],
            event_type=event['type'],
            processed=False
        )
        db.add(webhook_log)
        db.commit()

        # Handle different event types
        event_type = event['type']
        data = event['data']['object']

        if event_type == 'checkout.session.completed':
            # New subscription created
            user_id = data['metadata'].get('user_id')
            plan_id = data['metadata'].get('plan_id')

            if user_id and plan_id:
                subscription = Subscription(
                    user_id=user_id,
                    plan_id=int(plan_id),
                    stripe_customer_id=data['customer'],
                    stripe_subscription_id=data['subscription'],
                    status='active',
                    current_period_start=datetime.utcfromtimestamp(data['subscription_details']['current_period_start']) if 'subscription_details' in data else datetime.utcnow(),
                    current_period_end=datetime.utcfromtimestamp(data['subscription_details']['current_period_end']) if 'subscription_details' in data else datetime.utcnow()
                )
                db.add(subscription)

                # Audit log
                audit = AuditService(db)
                audit.log_action(
                    action="subscription_created",
                    user_id=user_id,
                    resource_type="subscription",
                    metadata={"plan_id": plan_id, "stripe_subscription_id": data['subscription']}
                )

        elif event_type == 'customer.subscription.updated':
            # Subscription updated (upgrade/downgrade/renewal)
            stripe_sub_id = data['id']
            subscription = db.query(Subscription).filter(
                Subscription.stripe_subscription_id == stripe_sub_id
            ).first()

            if subscription:
                subscription.status = data['status']
                subscription.current_period_start = datetime.utcfromtimestamp(data['current_period_start'])
                subscription.current_period_end = datetime.utcfromtimestamp(data['current_period_end'])
                subscription.cancel_at_period_end = data.get('cancel_at_period_end', False)

                audit = AuditService(db)
                audit.log_action(
                    action="subscription_updated",
                    user_id=subscription.user_id,
                    resource_type="subscription",
                    metadata={"status": data['status'], "stripe_subscription_id": stripe_sub_id}
                )

        elif event_type == 'customer.subscription.deleted':
            # Subscription cancelled/expired
            stripe_sub_id = data['id']
            subscription = db.query(Subscription).filter(
                Subscription.stripe_subscription_id == stripe_sub_id
            ).first()

            if subscription:
                subscription.status = 'cancelled'
                subscription.cancelled_at = datetime.utcnow()

                audit = AuditService(db)
                audit.log_action(
                    action="subscription_cancelled",
                    user_id=subscription.user_id,
                    resource_type="subscription",
                    metadata={"stripe_subscription_id": stripe_sub_id}
                )

        elif event_type == 'invoice.payment_failed':
            # Payment failed
            stripe_sub_id = data['subscription']
            subscription = db.query(Subscription).filter(
                Subscription.stripe_subscription_id == stripe_sub_id
            ).first()

            if subscription:
                subscription.status = 'past_due'

                audit = AuditService(db)
                audit.log_action(
                    action="payment_failed",
                    user_id=subscription.user_id,
                    resource_type="subscription",
                    metadata={"stripe_subscription_id": stripe_sub_id, "amount": data['amount_due']}
                )

        # Mark webhook as processed
        webhook_log.processed = True
        db.commit()

        return {"status": "success", "event_type": event_type}

    except Exception as e:
        db.rollback()
        # Log error but return 200 to prevent Stripe from retrying
        print(f"Webhook error: {str(e)}")
        return {"status": "error", "message": str(e)}
    finally:
        db.close()
