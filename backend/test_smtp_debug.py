import smtplib
import ssl
import os
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart

# Get credentials from .env
smtp_host = "smtp.hostinger.com"
smtp_port = 465
smtp_user = "robin.tomar@atgcflow.com"
smtp_password = "Rob@04051998"

print("=" * 60)
print("SMTP Connection Test - Detailed Diagnostics")
print("=" * 60)
print(f"\nSMTP Host: {smtp_host}")
print(f"SMTP Port: {smtp_port}")
print(f"SMTP User: {smtp_user}")
print(f"Password: {'*' * len(smtp_password)}")

print("\n" + "=" * 60)
print("Test 1: SSL Connection (Port 465)")
print("=" * 60)

try:
    context = ssl.create_default_context()
    print("Creating SSL connection...")
    server = smtplib.SMTP_SSL(smtp_host, smtp_port, context=context, timeout=10)
    print("✅ SSL connection established")
    
    print("Attempting login...")
    server.login(smtp_user, smtp_password)
    print("✅ Login successful!")
    
    server.quit()
    print("✅ Connection closed properly")
    print("\n🎉 SMTP Test PASSED!")
    
except smtplib.SMTPAuthenticationError as e:
    print(f"❌ Authentication failed: {e}")
    print("\nPossible fixes:")
    print("1. Verify password is correct")
    print("2. Check if 2FA is enabled - you may need app-specific password")
    print("3. Log into Hostinger email control panel and verify SMTP is enabled")
    print("4. Generate an app-specific password if available")
    
except Exception as e:
    print(f"❌ Connection failed: {type(e).__name__}: {e}")

print("\n" + "=" * 60)
print("Test 2: TLS Connection (Port 587)")
print("=" * 60)

try:
    context = ssl.create_default_context()
    print("Creating TLS connection...")
    server = smtplib.SMTP(smtp_host, 587, timeout=10)
    server.starttls(context=context)
    print("✅ TLS connection established")
    
    print("Attempting login...")
    server.login(smtp_user, smtp_password)
    print("✅ Login successful!")
    
    server.quit()
    print("✅ Connection closed properly")
    print("\n🎉 SMTP Test PASSED with TLS!")
    
except smtplib.SMTPAuthenticationError as e:
    print(f"❌ Authentication failed: {e}")
    
except Exception as e:
    print(f"❌ Connection failed: {type(e).__name__}: {e}")

print("\n" + "=" * 60)
print("Next Steps:")
print("=" * 60)
print("1. Log into Hostinger control panel (hpanel)")
print("2. Go to Emails > Email Accounts")
print("3. Check 'robin.tomar@atgcflow.com' settings")
print("4. Verify password or generate app-specific password")
print("5. Ensure SMTP access is enabled")
