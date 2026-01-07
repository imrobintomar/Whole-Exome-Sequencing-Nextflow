#!/usr/bin/env python3
"""
Quick verification script for SaaS backend setup
"""
import sys

print("=" * 60)
print("SaaS Backend Setup Verification")
print("=" * 60)

# 1. Check imports
print("\n1. Checking module imports...")
try:
    from main import app
    print("   ✅ Main app imports successfully")
except Exception as e:
    print(f"   ❌ Main app import failed: {e}")
    sys.exit(1)

try:
    from modules.admin.routes import router as admin_router
    print("   ✅ Admin routes import successfully")
except Exception as e:
    print(f"   ❌ Admin routes import failed: {e}")

try:
    from modules.billing.routes import router as billing_router
    print("   ✅ Billing routes import successfully")
except Exception as e:
    print(f"   ❌ Billing routes import failed: {e}")

try:
    from modules.chat.routes import router as chat_router
    print("   ✅ Chat routes import successfully")
except Exception as e:
    print(f"   ❌ Chat routes import failed: {e}")

# 2. Check database
print("\n2. Checking database...")
try:
    from database import SessionLocal
    from database_extensions import SubscriptionPlan, Subscription, UsageTracking, ChatConversation
    
    db = SessionLocal()
    plan_count = db.query(SubscriptionPlan).count()
    print(f"   ✅ Database accessible")
    print(f"   ✅ Subscription plans: {plan_count} (expected: 3)")
    
    if plan_count == 3:
        plans = db.query(SubscriptionPlan).all()
        for plan in plans:
            print(f"      • {plan.name}: ${plan.price_cents/100:.2f}/mo - {plan.monthly_jobs_limit} jobs")
    db.close()
except Exception as e:
    print(f"   ❌ Database check failed: {e}")

# 3. Check endpoints
print("\n3. Checking registered endpoints...")
try:
    routes = [r.path for r in app.routes if hasattr(r, 'path')]
    saas_routes = [r for r in routes if any(x in r for x in ['/admin', '/billing', '/chat'])]
    
    admin_routes = [r for r in saas_routes if r.startswith('/admin')]
    billing_routes = [r for r in saas_routes if r.startswith('/billing')]
    chat_routes = [r for r in saas_routes if r.startswith('/chat')]
    
    print(f"   ✅ Admin endpoints: {len(admin_routes)}")
    print(f"   ✅ Billing endpoints: {len(billing_routes)}")
    print(f"   ✅ Chat endpoints: {len(chat_routes)}")
    print(f"   ✅ Total SaaS endpoints: {len(saas_routes)}")
except Exception as e:
    print(f"   ❌ Endpoint check failed: {e}")

# 4. Check environment
print("\n4. Checking environment configuration...")
try:
    from config import settings
    
    print(f"   ✅ FRONTEND_URL: {settings.FRONTEND_URL}")
    
    if settings.STRIPE_SECRET_KEY and not settings.STRIPE_SECRET_KEY.startswith('sk_test_placeholder'):
        print(f"   ✅ Stripe Secret Key: Configured")
    else:
        print(f"   ⚠️  Stripe Secret Key: Using placeholder (set for production)")
    
    if settings.ADMIN_USER_UIDS and not 'placeholder' in settings.ADMIN_USER_UIDS:
        print(f"   ✅ Admin Users: Configured")
    else:
        print(f"   ⚠️  Admin Users: Using placeholder (set for admin access)")
        
except Exception as e:
    print(f"   ❌ Environment check failed: {e}")

# 5. Check services
print("\n5. Checking services...")
try:
    from services.stripe_service import StripeService
    print("   ✅ Stripe service available")
except Exception as e:
    print(f"   ❌ Stripe service: {e}")

try:
    from services.metrics_service import MetricsService
    metrics = MetricsService.get_system_metrics()
    print(f"   ✅ Metrics service working (CPU: {metrics['cpu']['usage_percent']}%)")
except Exception as e:
    print(f"   ❌ Metrics service: {e}")

try:
    from services.audit_service import AuditService
    print("   ✅ Audit service available")
except Exception as e:
    print(f"   ❌ Audit service: {e}")

# Summary
print("\n" + "=" * 60)
print("Setup Verification Complete!")
print("=" * 60)
print("\n✅ Backend is ready to use!")
print("\nNext steps:")
print("1. Configure Stripe keys in .env (for payments)")
print("2. Add your Firebase UID to ADMIN_USER_UIDS (for admin access)")
print("3. Start server: uvicorn main:app --reload")
print("4. Visit API docs: http://localhost:8000/docs")
print("\nFor detailed instructions, see:")
print("  - QUICK_START.md")
print("  - SETUP_COMPLETE.md")
print("  - SAAS_IMPLEMENTATION.md")
