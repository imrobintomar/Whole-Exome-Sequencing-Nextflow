#!/usr/bin/env python3
"""
Simple Email Test Script
Tests email sending with current configuration
"""

from dotenv import load_dotenv
load_dotenv()

from services.email_service import get_email_service

def main():
    print("=" * 70)
    print("TESTING EMAIL SYSTEM")
    print("=" * 70)

    email_service = get_email_service()

    # Test 1: Send to Hostinger email (same server)
    print("\n1. Testing email to robin.tomar@atgcflow.com...")
    print("-" * 70)
    success1 = email_service.send_custom_notification(
        user_email="robin.tomar@atgcflow.com",
        subject="Test Email - Same Server",
        message="This is a test email sent to the same Hostinger email server.",
        user_name="Robin"
    )

    if success1:
        print("Email sent successfully!")
        print("Check robin.tomar@atgcflow.com inbox")
    else:
        print(" Failed to send email")

    # Test 2: Send to admin email
    print("\n2. Testing email to admin@atgcflow.com...")
    print("-" * 70)
    success2 = email_service.send_custom_notification(
        user_email="admin@atgcflow.com",
        subject=" Test Email - Admin",
        message="This is a test email sent to the admin email address.",
        user_name="Admin"
    )

    if success2:
        print("Email sent successfully!")
        print("Check admin@atgcflow.com inbox")
    else:
        print("Failed to send email")

    print("\n" + "=" * 70)
    print("TEST COMPLETE")
    print("=" * 70)
    print("\nNote: If emails aren't in inbox, check:")
    print("  1. Spam/Junk folder")
    print("  2. Hostinger webmail at https://webmail.hostinger.com")
    print("  3. Email client settings")

if __name__ == "__main__":
    main()
