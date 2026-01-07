"""
Stripe Service
Wrapper for Stripe API operations
"""

import stripe
import os
from typing import Optional

# Initialize Stripe with API key from environment
stripe.api_key = os.getenv("STRIPE_SECRET_KEY")

class StripeService:
    def __init__(self):
        if not stripe.api_key:
            print("  Warning: STRIPE_SECRET_KEY not set. Billing features disabled.")

    def create_checkout_session(
        self,
        price_id: str,
        customer_email: str,
        success_url: str,
        cancel_url: str,
        metadata: dict = None
    ) -> Optional[str]:
        """
        Create Stripe Checkout session
        Returns: checkout session URL or None if Stripe not configured
        """
        if not stripe.api_key:
            return None

        try:
            session = stripe.checkout.Session.create(
                customer_email=customer_email,
                payment_method_types=['card'],
                line_items=[{
                    'price': price_id,
                    'quantity': 1,
                }],
                mode='subscription',
                success_url=success_url,
                cancel_url=cancel_url,
                metadata=metadata or {},
                allow_promotion_codes=True,
            )
            return session.url
        except stripe.error.StripeError as e:
            print(f"Stripe error: {e}")
            return None

    def create_customer_portal_session(
        self,
        customer_id: str,
        return_url: str
    ) -> Optional[str]:
        """
        Create Stripe Customer Portal session
        Returns: portal session URL or None if Stripe not configured
        """
        if not stripe.api_key:
            return None

        try:
            session = stripe.billing_portal.Session.create(
                customer=customer_id,
                return_url=return_url
            )
            return session.url
        except stripe.error.StripeError as e:
            print(f"Stripe error: {e}")
            return None

    def get_subscription(self, subscription_id: str):
        """Get subscription details from Stripe"""
        if not stripe.api_key:
            return None

        try:
            return stripe.Subscription.retrieve(subscription_id)
        except stripe.error.StripeError as e:
            print(f"Stripe error: {e}")
            return None

    def cancel_subscription(self, subscription_id: str, at_period_end: bool = True):
        """Cancel a subscription"""
        if not stripe.api_key:
            return None

        try:
            if at_period_end:
                return stripe.Subscription.modify(
                    subscription_id,
                    cancel_at_period_end=True
                )
            else:
                return stripe.Subscription.delete(subscription_id)
        except stripe.error.StripeError as e:
            print(f"Stripe error: {e}")
            return None

    def construct_webhook_event(self, payload: bytes, sig_header: str):
        """
        Construct and verify webhook event
        Raises stripe.error.SignatureVerificationError if invalid
        """
        webhook_secret = os.getenv("STRIPE_WEBHOOK_SECRET")
        if not webhook_secret:
            raise ValueError("STRIPE_WEBHOOK_SECRET not configured")

        return stripe.Webhook.construct_event(
            payload, sig_header, webhook_secret
        )
