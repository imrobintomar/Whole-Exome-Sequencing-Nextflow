"""
Email Notification Service - Handles all email notifications
Supports health alerts, payment notifications, and custom user emails
"""
import smtplib
import ssl
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from typing import List, Optional, Dict
from datetime import datetime
import os
from jinja2 import Template

class EmailService:
    """Service for sending email notifications"""

    def __init__(self):
        """Initialize email service with SMTP configuration"""
        self.smtp_host = os.getenv("SMTP_HOST", "smtp.hostinger.com")
        self.smtp_port = int(os.getenv("SMTP_PORT", "465"))
        self.smtp_user = os.getenv("SMTP_USER", "robin.tomar@atgcflow.com")
        self.smtp_password = os.getenv("SMTP_PASSWORD", "")
        self.admin_email = os.getenv("ADMIN_EMAIL", "admin@atgcflow.com")
        self.from_email = self.smtp_user
        self.from_name = os.getenv("EMAIL_FROM_NAME", "ATGC Flow")

    def _create_connection(self):
        """Create SMTP connection"""
        context = ssl.create_default_context()

        if self.smtp_port == 465:
            # SSL connection
            server = smtplib.SMTP_SSL(self.smtp_host, self.smtp_port, context=context)
        else:
            # TLS connection
            server = smtplib.SMTP(self.smtp_host, self.smtp_port)
            server.starttls(context=context)

        server.login(self.smtp_user, self.smtp_password)
        return server

    def _send_email(
        self,
        to_emails: List[str],
        subject: str,
        html_body: str,
        text_body: Optional[str] = None
    ) -> bool:
        """
        Send email using SMTP

        Args:
            to_emails: List of recipient email addresses
            subject: Email subject
            html_body: HTML email body
            text_body: Plain text email body (optional)

        Returns:
            bool: True if email sent successfully
        """
        try:
            message = MIMEMultipart("alternative")
            message["Subject"] = subject
            message["From"] = f"{self.from_name} <{self.from_email}>"
            message["To"] = ", ".join(to_emails)

            # Add plain text version
            if text_body:
                part1 = MIMEText(text_body, "plain")
                message.attach(part1)

            # Add HTML version
            part2 = MIMEText(html_body, "html")
            message.attach(part2)

            # Send email
            server = self._create_connection()
            server.sendmail(self.from_email, to_emails, message.as_string())
            server.quit()

            print(f"✅ Email sent successfully to {', '.join(to_emails)}")
            return True

        except Exception as e:
            print(f"❌ Failed to send email: {str(e)}")
            return False

    def send_health_alert(
        self,
        alert_type: str,
        severity: str,
        message: str,
        value: float,
        threshold: float,
        metrics: Optional[Dict] = None
    ) -> bool:
        """
        Send health alert email to admin

        Args:
            alert_type: Type of alert (cpu, memory, disk)
            severity: Alert severity (warning, critical)
            message: Alert message
            value: Current value
            threshold: Threshold value
            metrics: Additional metrics data

        Returns:
            bool: True if email sent successfully
        """
        subject = f"🚨 {severity.upper()} System Alert: {alert_type.upper()}"

        # Create HTML email body
        html_body = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <style>
                body {{ font-family: Arial, sans-serif; line-height: 1.6; color: #333; }}
                .container {{ max-width: 600px; margin: 0 auto; padding: 20px; }}
                .header {{ background-color: {'#dc2626' if severity == 'critical' else '#f59e0b'};
                           color: white; padding: 20px; border-radius: 8px 8px 0 0; }}
                .content {{ background-color: #f9fafb; padding: 20px; border-radius: 0 0 8px 8px; }}
                .alert-box {{ background-color: white; padding: 15px; margin: 15px 0;
                             border-left: 4px solid {'#dc2626' if severity == 'critical' else '#f59e0b'}; }}
                .metric {{ display: flex; justify-content: space-between; padding: 8px 0;
                          border-bottom: 1px solid #e5e7eb; }}
                .metric:last-child {{ border-bottom: none; }}
                .label {{ font-weight: bold; color: #6b7280; }}
                .value {{ color: #111827; }}
                .footer {{ margin-top: 20px; padding-top: 20px; border-top: 1px solid #e5e7eb;
                          font-size: 12px; color: #6b7280; text-align: center; }}
            </style>
        </head>
        <body>
            <div class="container">
                <div class="header">
                    <h2 style="margin: 0;">System Health Alert</h2>
                    <p style="margin: 5px 0 0 0; opacity: 0.9;">ATGC Flow - Whole Exome Sequencing Platform</p>
                </div>
                <div class="content">
                    <div class="alert-box">
                        <h3 style="margin-top: 0; color: {'#dc2626' if severity == 'critical' else '#f59e0b'};">
                            {severity.upper()} Alert
                        </h3>
                        <p style="font-size: 16px; margin: 10px 0;">{message}</p>
                        <div style="margin-top: 15px;">
                            <div class="metric">
                                <span class="label">Alert Type:</span>
                                <span class="value">{alert_type.upper()}</span>
                            </div>
                            <div class="metric">
                                <span class="label">Current Value:</span>
                                <span class="value">{value:.1f}%</span>
                            </div>
                            <div class="metric">
                                <span class="label">Threshold:</span>
                                <span class="value">{threshold:.1f}%</span>
                            </div>
                            <div class="metric">
                                <span class="label">Timestamp:</span>
                                <span class="value">{datetime.now().strftime('%Y-%m-%d %H:%M:%S UTC')}</span>
                            </div>
                        </div>
                    </div>
        """

        # Add system metrics if provided
        if metrics:
            html_body += """
                    <div style="margin-top: 20px;">
                        <h4 style="margin-bottom: 10px; color: #374151;">System Metrics</h4>
                        <div style="background-color: white; padding: 15px; border-radius: 8px;">
            """
            for key, val in metrics.items():
                html_body += f"""
                            <div class="metric">
                                <span class="label">{key.upper()}:</span>
                                <span class="value">{val:.1f}%</span>
                            </div>
                """
            html_body += """
                        </div>
                    </div>
            """

        html_body += """
                    <p style="margin-top: 20px; color: #6b7280;">
                        Please check the admin dashboard for more details and take appropriate action.
                    </p>
                </div>
                <div class="footer">
                    <p>This is an automated alert from ATGC Flow Health Monitoring System</p>
                    <p>&copy; 2026 ATGC Flow. All rights reserved.</p>
                </div>
            </div>
        </body>
        </html>
        """

        # Plain text version
        text_body = f"""
        System Health Alert - {severity.upper()}

        {message}

        Alert Type: {alert_type.upper()}
        Current Value: {value:.1f}%
        Threshold: {threshold:.1f}%
        Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S UTC')}

        Please check the admin dashboard for more details.

        ---
        ATGC Flow - Whole Exome Sequencing Platform
        """

        return self._send_email([self.admin_email], subject, html_body, text_body)

    def send_payment_due_alert(
        self,
        user_email: str,
        user_name: str,
        amount: float,
        due_date: str,
        invoice_url: Optional[str] = None
    ) -> bool:
        """
        Send payment due notification to user

        Args:
            user_email: User email address
            user_name: User name
            amount: Payment amount
            due_date: Payment due date
            invoice_url: Optional invoice URL

        Returns:
            bool: True if email sent successfully
        """
        subject = f"Payment Due Reminder - ATGC Flow"

        html_body = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <style>
                body {{ font-family: Arial, sans-serif; line-height: 1.6; color: #333; }}
                .container {{ max-width: 600px; margin: 0 auto; padding: 20px; }}
                .header {{ background-color: #3b82f6; color: white; padding: 20px; border-radius: 8px 8px 0 0; }}
                .content {{ background-color: #f9fafb; padding: 20px; border-radius: 0 0 8px 8px; }}
                .info-box {{ background-color: white; padding: 20px; margin: 15px 0; border-radius: 8px; }}
                .amount {{ font-size: 32px; font-weight: bold; color: #3b82f6; text-align: center; margin: 20px 0; }}
                .button {{ display: inline-block; background-color: #3b82f6; color: white;
                          padding: 12px 30px; text-decoration: none; border-radius: 6px;
                          font-weight: bold; margin: 20px 0; }}
                .footer {{ margin-top: 20px; padding-top: 20px; border-top: 1px solid #e5e7eb;
                          font-size: 12px; color: #6b7280; text-align: center; }}
            </style>
        </head>
        <body>
            <div class="container">
                <div class="header">
                    <h2 style="margin: 0;">Payment Reminder</h2>
                    <p style="margin: 5px 0 0 0; opacity: 0.9;">ATGC Flow</p>
                </div>
                <div class="content">
                    <p>Hi {user_name},</p>
                    <p>This is a friendly reminder that your payment is due soon.</p>

                    <div class="info-box">
                        <div class="amount">${amount:.2f}</div>
                        <p style="text-align: center; color: #6b7280;">Due Date: <strong>{due_date}</strong></p>
                    </div>

                    <p>Please ensure your payment is processed by the due date to avoid any service interruption.</p>
        """

        if invoice_url:
            html_body += f"""
                    <div style="text-align: center;">
                        <a href="{invoice_url}" class="button">View Invoice</a>
                    </div>
            """

        html_body += """
                    <p style="margin-top: 20px;">If you have any questions or concerns, please don't hesitate to contact our support team.</p>

                    <p>Thank you for using ATGC Flow!</p>
                </div>
                <div class="footer">
                    <p>ATGC Flow - Whole Exome Sequencing Platform</p>
                    <p>&copy; 2026 ATGC Flow. All rights reserved.</p>
                </div>
            </div>
        </body>
        </html>
        """

        text_body = f"""
        Payment Reminder

        Hi {user_name},

        This is a friendly reminder that your payment is due soon.

        Amount Due: ${amount:.2f}
        Due Date: {due_date}

        Please ensure your payment is processed by the due date to avoid any service interruption.

        {"View Invoice: " + invoice_url if invoice_url else ""}

        Thank you for using ATGC Flow!

        ---
        ATGC Flow - Whole Exome Sequencing Platform
        """

        return self._send_email([user_email], subject, html_body, text_body)

    def send_subscription_expiry_alert(
        self,
        user_email: str,
        user_name: str,
        expiry_date: str,
        plan_name: str
    ) -> bool:
        """
        Send subscription expiry notification to user

        Args:
            user_email: User email address
            user_name: User name
            expiry_date: Subscription expiry date
            plan_name: Current plan name

        Returns:
            bool: True if email sent successfully
        """
        subject = f"Subscription Expiring Soon - ATGC Flow"

        html_body = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <style>
                body {{ font-family: Arial, sans-serif; line-height: 1.6; color: #333; }}
                .container {{ max-width: 600px; margin: 0 auto; padding: 20px; }}
                .header {{ background-color: #f59e0b; color: white; padding: 20px; border-radius: 8px 8px 0 0; }}
                .content {{ background-color: #f9fafb; padding: 20px; border-radius: 0 0 8px 8px; }}
                .info-box {{ background-color: white; padding: 20px; margin: 15px 0; border-radius: 8px;
                            border-left: 4px solid #f59e0b; }}
                .button {{ display: inline-block; background-color: #f59e0b; color: white;
                          padding: 12px 30px; text-decoration: none; border-radius: 6px;
                          font-weight: bold; margin: 20px 0; }}
                .footer {{ margin-top: 20px; padding-top: 20px; border-top: 1px solid #e5e7eb;
                          font-size: 12px; color: #6b7280; text-align: center; }}
            </style>
        </head>
        <body>
            <div class="container">
                <div class="header">
                    <h2 style="margin: 0;">Subscription Expiring Soon</h2>
                    <p style="margin: 5px 0 0 0; opacity: 0.9;">ATGC Flow</p>
                </div>
                <div class="content">
                    <p>Hi {user_name},</p>
                    <p>Your <strong>{plan_name}</strong> subscription is expiring soon.</p>

                    <div class="info-box">
                        <p style="margin: 0; font-size: 18px;">
                            <strong>Expiry Date:</strong> {expiry_date}
                        </p>
                    </div>

                    <p>To continue enjoying uninterrupted access to our services, please renew your subscription.</p>

                    <div style="text-align: center;">
                        <a href="{os.getenv('FRONTEND_URL', 'http://localhost:3000')}/settings" class="button">
                            Renew Subscription
                        </a>
                    </div>

                    <p style="margin-top: 20px;">If you have any questions, please contact our support team.</p>

                    <p>Thank you for using ATGC Flow!</p>
                </div>
                <div class="footer">
                    <p>ATGC Flow - Whole Exome Sequencing Platform</p>
                    <p>&copy; 2026 ATGC Flow. All rights reserved.</p>
                </div>
            </div>
        </body>
        </html>
        """

        text_body = f"""
        Subscription Expiring Soon

        Hi {user_name},

        Your {plan_name} subscription is expiring soon.

        Expiry Date: {expiry_date}

        To continue enjoying uninterrupted access to our services, please renew your subscription.

        Visit: {os.getenv('FRONTEND_URL', 'http://localhost:3000')}/settings

        Thank you for using ATGC Flow!

        ---
        ATGC Flow - Whole Exome Sequencing Platform
        """

        return self._send_email([user_email], subject, html_body, text_body)

    def send_custom_notification(
        self,
        user_email: str,
        subject: str,
        message: str,
        user_name: Optional[str] = None
    ) -> bool:
        """
        Send custom notification email

        Args:
            user_email: User email address
            subject: Email subject
            message: Email message
            user_name: Optional user name

        Returns:
            bool: True if email sent successfully
        """
        greeting = f"Hi {user_name}," if user_name else "Hello,"

        html_body = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <style>
                body {{ font-family: Arial, sans-serif; line-height: 1.6; color: #333; }}
                .container {{ max-width: 600px; margin: 0 auto; padding: 20px; }}
                .header {{ background-color: #6366f1; color: white; padding: 20px; border-radius: 8px 8px 0 0; }}
                .content {{ background-color: #f9fafb; padding: 20px; border-radius: 0 0 8px 8px; }}
                .message-box {{ background-color: white; padding: 20px; margin: 15px 0; border-radius: 8px; }}
                .footer {{ margin-top: 20px; padding-top: 20px; border-top: 1px solid #e5e7eb;
                          font-size: 12px; color: #6b7280; text-align: center; }}
            </style>
        </head>
        <body>
            <div class="container">
                <div class="header">
                    <h2 style="margin: 0;">ATGC Flow Notification</h2>
                </div>
                <div class="content">
                    <p>{greeting}</p>

                    <div class="message-box">
                        <p style="white-space: pre-wrap;">{message}</p>
                    </div>

                    <p>If you have any questions, please contact our support team.</p>

                    <p>Best regards,<br>ATGC Flow Team</p>
                </div>
                <div class="footer">
                    <p>ATGC Flow - Whole Exome Sequencing Platform</p>
                    <p>&copy; 2026 ATGC Flow. All rights reserved.</p>
                </div>
            </div>
        </body>
        </html>
        """

        text_body = f"""
        {greeting}

        {message}

        If you have any questions, please contact our support team.

        Best regards,
        ATGC Flow Team

        ---
        ATGC Flow - Whole Exome Sequencing Platform
        """

        return self._send_email([user_email], subject, html_body, text_body)

    def test_connection(self) -> Dict:
        """
        Test SMTP connection and configuration

        Returns:
            Dict with connection status and details
        """
        try:
            server = self._create_connection()
            server.quit()
            return {
                "status": "success",
                "message": "SMTP connection successful",
                "smtp_host": self.smtp_host,
                "smtp_port": self.smtp_port,
                "smtp_user": self.smtp_user,
                "admin_email": self.admin_email
            }
        except Exception as e:
            return {
                "status": "error",
                "message": f"SMTP connection failed: {str(e)}",
                "smtp_host": self.smtp_host,
                "smtp_port": self.smtp_port
            }


# Singleton instance
_email_service = None

def get_email_service() -> EmailService:
    """Get or create email service singleton"""
    global _email_service
    if _email_service is None:
        _email_service = EmailService()
    return _email_service
