"""
Database Extensions for SaaS Features
ADDITIVE ONLY - Does not modify existing tables
"""

from sqlalchemy import Column, Integer, String, Boolean, DateTime, Text, ForeignKey, Index
from sqlalchemy.orm import relationship
from datetime import datetime
from database import Base

# ============================================================================
# SUBSCRIPTION MANAGEMENT
# ============================================================================

class SubscriptionPlan(Base):
    __tablename__ = 'subscription_plans'

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(50), nullable=False, unique=True)
    stripe_price_id = Column(String(100), nullable=True)
    monthly_jobs_limit = Column(Integer, nullable=False)
    chat_support = Column(Boolean, default=False)
    price_cents = Column(Integer, nullable=False)
    features_json = Column(Text, nullable=True)
    active = Column(Boolean, default=True)
    created_at = Column(DateTime, default=datetime.utcnow)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)


class Subscription(Base):
    __tablename__ = 'subscriptions'

    id = Column(Integer, primary_key=True, autoincrement=True)
    user_id = Column(String(128), nullable=False)
    plan_id = Column(Integer, ForeignKey('subscription_plans.id'), nullable=False)
    stripe_customer_id = Column(String(100), nullable=True)
    stripe_subscription_id = Column(String(100), nullable=True)
    status = Column(String(20), nullable=False, default='active')
    current_period_start = Column(DateTime, nullable=True)
    current_period_end = Column(DateTime, nullable=True)
    cancel_at_period_end = Column(Boolean, default=False)
    created_at = Column(DateTime, default=datetime.utcnow)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)

    # Relationship
    plan = relationship("SubscriptionPlan", backref="subscriptions")


class UsageTracking(Base):
    __tablename__ = 'usage_tracking'

    id = Column(Integer, primary_key=True, autoincrement=True)
    user_id = Column(String(128), nullable=False)
    month = Column(Integer, nullable=False)  # YYYYMM format
    jobs_executed = Column(Integer, default=0)
    jobs_limit = Column(Integer, nullable=False)
    created_at = Column(DateTime, default=datetime.utcnow)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)

    __table_args__ = (
        Index('idx_usage_user_month', 'user_id', 'month', unique=True),
    )


# ============================================================================
# CHAT SYSTEM
# ============================================================================

class ChatConversation(Base):
    __tablename__ = 'chat_conversations'

    id = Column(Integer, primary_key=True, autoincrement=True)
    user_id = Column(String(128), nullable=False)
    job_id = Column(String(36), nullable=True)
    subject = Column(String(200), nullable=True)
    status = Column(String(20), default='open')  # open, resolved, closed
    last_message_at = Column(DateTime, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)

    # Relationship
    messages = relationship("ChatMessage", back_populates="conversation", cascade="all, delete-orphan")


class ChatMessage(Base):
    __tablename__ = 'chat_messages'

    id = Column(Integer, primary_key=True, autoincrement=True)
    conversation_id = Column(Integer, ForeignKey('chat_conversations.id'), nullable=False)
    sender_id = Column(String(128), nullable=False)
    sender_role = Column(String(10), nullable=False)  # user, admin
    message_text = Column(Text, nullable=False)
    read_by_admin = Column(Boolean, default=False)
    read_by_user = Column(Boolean, default=False)
    created_at = Column(DateTime, default=datetime.utcnow)

    # Relationship
    conversation = relationship("ChatConversation", back_populates="messages")


# ============================================================================
# ADMIN & AUDIT
# ============================================================================

class AdminUser(Base):
    __tablename__ = 'admin_users'

    id = Column(Integer, primary_key=True, autoincrement=True)
    firebase_uid = Column(String(128), unique=True, nullable=False)
    email = Column(String(255), nullable=False)
    is_super_admin = Column(Boolean, default=False)
    created_at = Column(DateTime, default=datetime.utcnow)
    last_login = Column(DateTime, nullable=True)


class AuditLog(Base):
    __tablename__ = 'audit_logs'

    id = Column(Integer, primary_key=True, autoincrement=True)
    user_id = Column(String(128), nullable=True)
    action = Column(String(100), nullable=False)
    resource_type = Column(String(50), nullable=True)
    resource_id = Column(String(100), nullable=True)
    ip_address = Column(String(45), nullable=True)
    user_agent = Column(Text, nullable=True)
    metadata_json = Column(Text, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow)

    __table_args__ = (
        Index('idx_audit_user', 'user_id'),
        Index('idx_audit_action', 'action'),
    )


class WebhookEvent(Base):
    __tablename__ = 'webhook_events'

    id = Column(Integer, primary_key=True, autoincrement=True)
    stripe_event_id = Column(String(100), unique=True, nullable=False)
    event_type = Column(String(100), nullable=False)
    processed = Column(Boolean, default=False)
    payload_json = Column(Text, nullable=False)
    processing_error = Column(Text, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow)
    processed_at = Column(DateTime, nullable=True)

    __table_args__ = (
        Index('idx_webhook_processed', 'processed'),
    )


# ============================================================================
# INDEXES FOR PERFORMANCE
# ============================================================================

Index('idx_subscriptions_user', Subscription.user_id)
Index('idx_subscriptions_stripe', Subscription.stripe_subscription_id)
Index('idx_conversations_user', ChatConversation.user_id)
Index('idx_conversations_status', ChatConversation.status)
Index('idx_messages_conversation', ChatMessage.conversation_id)
Index('idx_messages_unread_admin', ChatMessage.read_by_admin, ChatMessage.sender_role)
