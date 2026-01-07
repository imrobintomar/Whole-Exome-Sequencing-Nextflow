"""
Initialize SaaS Database Tables
Run this once to create new tables
"""

from database import engine
from database_extensions import Base, SubscriptionPlan
from sqlalchemy.orm import Session
import json

def init_saas_tables():
    """Create all SaaS extension tables"""
    print("Creating SaaS extension tables...")
    Base.metadata.create_all(bind=engine)
    print("✅ Tables created successfully")

def seed_subscription_plans():
    """Seed initial subscription plans"""
    print("Seeding subscription plans...")

    from database import SessionLocal
    db = SessionLocal()

    try:
        # Check if plans already exist
        existing = db.query(SubscriptionPlan).count()
        if existing > 0:
            print("⚠️  Plans already exist, skipping seed")
            return

        plans = [
            SubscriptionPlan(
                name="Free",
                stripe_price_id=None,
                monthly_jobs_limit=2,
                chat_support=False,
                price_cents=0,
                features_json=json.dumps([
                    "2 WES jobs per month",
                    "Email support",
                    "48-hour turnaround"
                ]),
                active=True
            ),
            SubscriptionPlan(
                name="Basic",
                stripe_price_id=None,  # Will be set via ENV
                monthly_jobs_limit=10,
                chat_support=True,
                price_cents=2900,
                features_json=json.dumps([
                    "10 WES jobs per month",
                    "Live chat support",
                    "24-hour turnaround",
                    "Priority processing"
                ]),
                active=True
            ),
            SubscriptionPlan(
                name="Pro",
                stripe_price_id=None,  # Will be set via ENV
                monthly_jobs_limit=50,
                chat_support=True,
                price_cents=9900,
                features_json=json.dumps([
                    "50 WES jobs per month",
                    "Live chat support",
                    "12-hour turnaround",
                    "Priority processing",
                    "Dedicated support",
                    "API access"
                ]),
                active=True
            )
        ]

        for plan in plans:
            db.add(plan)

        db.commit()
        print("✅ Subscription plans seeded")

        # Display plans
        for plan in plans:
            print(f"   • {plan.name}: ${plan.price_cents/100:.2f}/month - {plan.monthly_jobs_limit} jobs")

    except Exception as e:
        print(f"❌ Error seeding plans: {e}")
        db.rollback()
    finally:
        db.close()

if __name__ == "__main__":
    print("=" * 60)
    print("WES Platform - SaaS Database Initialization")
    print("=" * 60)

    init_saas_tables()
    seed_subscription_plans()

    print("\n" + "=" * 60)
    print("✅ SaaS database initialization complete!")
    print("=" * 60)
    print("\nNext steps:")
    print("1. Set required environment variables (see .env.saas.example)")
    print("2. Update Stripe Price IDs in subscription_plans table")
    print("3. Add admin users to admin_users table")
    print("4. Restart backend server")
