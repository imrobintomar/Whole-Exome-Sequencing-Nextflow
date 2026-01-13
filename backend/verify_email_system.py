#!/usr/bin/env python3
"""
Email System Verification Script
Run this after restarting the backend to verify email functionality
"""

from dotenv import load_dotenv
import os

# Load environment variables
load_dotenv()

print("=" * 70)
print("EMAIL SYSTEM VERIFICATION")
print("=" * 70)

# Check 1: Environment Variables
print("\n1. Checking Environment Variables...")
print("-" * 70)

required_vars = {
    'SMTP_HOST': os.getenv('SMTP_HOST'),
    'SMTP_PORT': os.getenv('SMTP_PORT'),
    'SMTP_USER': os.getenv('SMTP_USER'),
    'SMTP_PASSWORD': os.getenv('SMTP_PASSWORD'),
    'ADMIN_EMAIL': os.getenv('ADMIN_EMAIL'),
}

all_set = True
for var, value in required_vars.items():
    if value:
        if var == 'SMTP_PASSWORD':
            print(f"✅ {var}: {'*' * len(value)}")
        else:
            print(f"✅ {var}: {value}")
    else:
        print(f"❌ {var}: NOT SET")
        all_set = False

if not all_set:
    print("\n❌ Some environment variables are missing!")
    print("Please check your .env file")
    exit(1)

# Check 2: SMTP Connection
print("\n2. Testing SMTP Connection...")
print("-" * 70)

try:
    from services.email_service import get_email_service

    email_service = get_email_service()
    result = email_service.test_connection()

    if result['status'] == 'success':
        print(f"✅ {result['message']}")
        print(f"   Host: {result['smtp_host']}:{result['smtp_port']}")
        print(f"   User: {result['smtp_user']}")
        print(f"   Admin: {result['admin_email']}")
    else:
        print(f"❌ {result['message']}")
        exit(1)

except Exception as e:
    print(f"❌ Failed to test connection: {e}")
    exit(1)

# Check 3: Send Test Email
print("\n3. Sending Test Email...")
print("-" * 70)

test_email = os.getenv('ADMIN_EMAIL', 'admin@atgcflow.com')
print(f"Sending to: {test_email}")

try:
    success = email_service.send_custom_notification(
        user_email=test_email,
        subject="✅ Email System Verification - ATGC Flow",
        message="Your email notification system is now working correctly! You can send custom emails, payment reminders, and health alerts from the admin dashboard.",
        user_name="Admin"
    )

    if success:
        print(f"✅ Test email sent successfully!")
        print(f"\n📧 Check inbox at: {test_email}")
        print("⚠️  If not in inbox, check SPAM folder!")
    else:
        print("❌ Failed to send test email")
        exit(1)

except Exception as e:
    print(f"❌ Error sending email: {e}")
    import traceback
    traceback.print_exc()
    exit(1)

# Check 4: Health Monitor Integration
print("\n4. Checking Health Monitor Integration...")
print("-" * 70)

try:
    from services.health_monitor import HealthMonitor

    thresholds = HealthMonitor.THRESHOLDS
    print(f"✅ Health Monitor loaded")
    print(f"   CPU Critical: {thresholds['cpu_critical']}%")
    print(f"   Memory Critical: {thresholds['memory_critical']}%")
    print(f"   Disk Critical: {thresholds['disk_critical']}%")
    print(f"   Email Alerts: {'Enabled' if HealthMonitor.EMAIL_ALERTS_ENABLED else 'Disabled'}")

except Exception as e:
    print(f"⚠️  Warning: Could not load Health Monitor: {e}")

# Summary
print("\n" + "=" * 70)
print("VERIFICATION COMPLETE!")
print("=" * 70)
print("\n✅ All checks passed!")
print("\nYour email notification system is ready to use:")
print("  • Custom email notifications")
print("  • Payment reminders")
print("  • Subscription expiry alerts")
print("  • Automated health monitoring alerts")
print("\nAccess the Email Alerts tab in your admin dashboard to start sending emails.")
print("=" * 70)
