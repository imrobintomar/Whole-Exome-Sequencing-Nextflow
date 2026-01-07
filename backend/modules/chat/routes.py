"""
Chat API Routes
User-Admin chat support system
"""

from fastapi import APIRouter, Depends, HTTPException
from firebase_auth import get_current_user
from middleware.admin_guard import require_admin
from database import SessionLocal, User
from database_extensions import ChatConversation, ChatMessage, Subscription, SubscriptionPlan
from services.audit_service import AuditService
from datetime import datetime
from pydantic import BaseModel

router = APIRouter(prefix="/chat", tags=["chat"])


class CreateConversationRequest(BaseModel):
    subject: str
    initial_message: str


class SendMessageRequest(BaseModel):
    message: str


# ============================================================================
# USER ENDPOINTS
# ============================================================================

@router.post("/conversations")
async def create_conversation(
    data: CreateConversationRequest,
    current_user=Depends(get_current_user)
):
    """
    Create new chat conversation (user)
    Requires active subscription with chat support
    """
    db = SessionLocal()
    try:
        # Check if user has chat support enabled
        subscription = db.query(Subscription).filter(
            Subscription.user_id == current_user.firebase_uid,
            Subscription.status.in_(['active', 'trialing'])
        ).first()

        if subscription:
            plan = db.query(SubscriptionPlan).filter(
                SubscriptionPlan.id == subscription.plan_id
            ).first()

            if not plan or not plan.chat_support:
                raise HTTPException(
                    status_code=403,
                    detail="Chat support not available on your plan. Please upgrade."
                )
        else:
            raise HTTPException(
                status_code=403,
                detail="Chat support requires an active subscription. Please upgrade."
            )

        # Create conversation
        conversation = ChatConversation(
            user_id=current_user.firebase_uid,
            subject=data.subject,
            status='open'
        )
        db.add(conversation)
        db.flush()

        # Add initial message
        message = ChatMessage(
            conversation_id=conversation.id,
            sender_id=current_user.firebase_uid,
            sender_type='user',
            message=data.initial_message
        )
        db.add(message)

        # Update conversation last message time
        conversation.last_message_at = datetime.utcnow()
        db.commit()

        # Audit log
        audit = AuditService(db)
        audit.log_action(
            action="chat_conversation_created",
            user_id=current_user.firebase_uid,
            resource_type="chat",
            resource_id=str(conversation.id),
            metadata={"subject": data.subject}
        )

        return {
            "conversation_id": conversation.id,
            "subject": conversation.subject,
            "status": conversation.status,
            "created_at": conversation.created_at
        }
    finally:
        db.close()


@router.get("/conversations")
async def get_my_conversations(current_user=Depends(get_current_user)):
    """
    Get all conversations for current user
    """
    db = SessionLocal()
    try:
        conversations = db.query(ChatConversation).filter(
            ChatConversation.user_id == current_user.firebase_uid
        ).order_by(ChatConversation.last_message_at.desc()).all()

        return {
            "conversations": [
                {
                    "id": conv.id,
                    "subject": conv.subject,
                    "status": conv.status,
                    "created_at": conv.created_at,
                    "last_message_at": conv.last_message_at,
                    "admin_id": conv.admin_id
                }
                for conv in conversations
            ]
        }
    finally:
        db.close()


@router.get("/conversations/{conversation_id}/messages")
async def get_conversation_messages(
    conversation_id: int,
    current_user=Depends(get_current_user)
):
    """
    Get all messages in a conversation
    """
    db = SessionLocal()
    try:
        # Verify user owns this conversation
        conversation = db.query(ChatConversation).filter(
            ChatConversation.id == conversation_id,
            ChatConversation.user_id == current_user.firebase_uid
        ).first()

        if not conversation:
            raise HTTPException(status_code=404, detail="Conversation not found")

        messages = db.query(ChatMessage).filter(
            ChatMessage.conversation_id == conversation_id
        ).order_by(ChatMessage.created_at.asc()).all()

        return {
            "conversation": {
                "id": conversation.id,
                "subject": conversation.subject,
                "status": conversation.status
            },
            "messages": [
                {
                    "id": msg.id,
                    "sender_type": msg.sender_type,
                    "message": msg.message,
                    "created_at": msg.created_at
                }
                for msg in messages
            ]
        }
    finally:
        db.close()


@router.post("/conversations/{conversation_id}/messages")
async def send_message(
    conversation_id: int,
    data: SendMessageRequest,
    current_user=Depends(get_current_user)
):
    """
    Send message in conversation (user)
    """
    db = SessionLocal()
    try:
        # Verify user owns this conversation
        conversation = db.query(ChatConversation).filter(
            ChatConversation.id == conversation_id,
            ChatConversation.user_id == current_user.firebase_uid
        ).first()

        if not conversation:
            raise HTTPException(status_code=404, detail="Conversation not found")

        if conversation.status == 'closed':
            raise HTTPException(
                status_code=400,
                detail="Cannot send messages to closed conversations"
            )

        # Create message
        message = ChatMessage(
            conversation_id=conversation_id,
            sender_id=current_user.firebase_uid,
            sender_type='user',
            message=data.message
        )
        db.add(message)

        # Update conversation
        conversation.last_message_at = datetime.utcnow()
        conversation.status = 'open'  # Reopen if was resolved
        db.commit()

        return {
            "message_id": message.id,
            "created_at": message.created_at
        }
    finally:
        db.close()


@router.patch("/conversations/{conversation_id}/close")
async def close_conversation(
    conversation_id: int,
    current_user=Depends(get_current_user)
):
    """
    Close/resolve conversation (user)
    """
    db = SessionLocal()
    try:
        conversation = db.query(ChatConversation).filter(
            ChatConversation.id == conversation_id,
            ChatConversation.user_id == current_user.firebase_uid
        ).first()

        if not conversation:
            raise HTTPException(status_code=404, detail="Conversation not found")

        conversation.status = 'closed'
        db.commit()

        # Audit log
        audit = AuditService(db)
        audit.log_action(
            action="chat_conversation_closed",
            user_id=current_user.firebase_uid,
            resource_type="chat",
            resource_id=str(conversation_id)
        )

        return {"status": "closed"}
    finally:
        db.close()


# ============================================================================
# ADMIN ENDPOINTS
# ============================================================================

@router.get("/admin/conversations")
async def get_all_conversations(
    status: str = None,
    limit: int = 50,
    offset: int = 0,
    admin=Depends(require_admin)
):
    """
    Get all chat conversations (admin)
    """
    db = SessionLocal()
    try:
        query = db.query(ChatConversation)

        if status:
            query = query.filter(ChatConversation.status == status)

        total = query.count()
        conversations = query.order_by(
            ChatConversation.last_message_at.desc()
        ).limit(limit).offset(offset).all()

        # Get user emails
        conv_data = []
        for conv in conversations:
            user = db.query(User).filter(User.firebase_uid == conv.user_id).first()
            conv_data.append({
                "id": conv.id,
                "subject": conv.subject,
                "status": conv.status,
                "user_id": conv.user_id,
                "user_email": user.email if user else None,
                "admin_id": conv.admin_id,
                "created_at": conv.created_at,
                "last_message_at": conv.last_message_at
            })

        return {
            "conversations": conv_data,
            "total": total,
            "limit": limit,
            "offset": offset
        }
    finally:
        db.close()


@router.get("/admin/conversations/{conversation_id}/messages")
async def admin_get_conversation_messages(
    conversation_id: int,
    admin=Depends(require_admin)
):
    """
    Get conversation messages (admin)
    """
    db = SessionLocal()
    try:
        conversation = db.query(ChatConversation).filter(
            ChatConversation.id == conversation_id
        ).first()

        if not conversation:
            raise HTTPException(status_code=404, detail="Conversation not found")

        # Get user info
        user = db.query(User).filter(User.firebase_uid == conversation.user_id).first()

        messages = db.query(ChatMessage).filter(
            ChatMessage.conversation_id == conversation_id
        ).order_by(ChatMessage.created_at.asc()).all()

        return {
            "conversation": {
                "id": conversation.id,
                "subject": conversation.subject,
                "status": conversation.status,
                "user_id": conversation.user_id,
                "user_email": user.email if user else None
            },
            "messages": [
                {
                    "id": msg.id,
                    "sender_type": msg.sender_type,
                    "sender_id": msg.sender_id,
                    "message": msg.message,
                    "created_at": msg.created_at
                }
                for msg in messages
            ]
        }
    finally:
        db.close()


@router.post("/admin/conversations/{conversation_id}/messages")
async def admin_send_message(
    conversation_id: int,
    data: SendMessageRequest,
    admin=Depends(require_admin)
):
    """
    Send message as admin
    """
    db = SessionLocal()
    try:
        conversation = db.query(ChatConversation).filter(
            ChatConversation.id == conversation_id
        ).first()

        if not conversation:
            raise HTTPException(status_code=404, detail="Conversation not found")

        # Assign admin to conversation if not assigned
        if not conversation.admin_id:
            conversation.admin_id = admin.firebase_uid

        # Create message
        message = ChatMessage(
            conversation_id=conversation_id,
            sender_id=admin.firebase_uid,
            sender_type='admin',
            message=data.message
        )
        db.add(message)

        # Update conversation
        conversation.last_message_at = datetime.utcnow()
        db.commit()

        # Audit log
        audit = AuditService(db)
        audit.log_action(
            action="admin_chat_message_sent",
            user_id=admin.firebase_uid,
            resource_type="chat",
            resource_id=str(conversation_id),
            metadata={"user_id": conversation.user_id}
        )

        return {
            "message_id": message.id,
            "created_at": message.created_at
        }
    finally:
        db.close()


@router.patch("/admin/conversations/{conversation_id}/status")
async def admin_update_conversation_status(
    conversation_id: int,
    status: str,
    admin=Depends(require_admin)
):
    """
    Update conversation status (admin)
    Valid statuses: open, resolved, closed
    """
    db = SessionLocal()
    try:
        if status not in ['open', 'resolved', 'closed']:
            raise HTTPException(
                status_code=400,
                detail="Invalid status. Must be: open, resolved, or closed"
            )

        conversation = db.query(ChatConversation).filter(
            ChatConversation.id == conversation_id
        ).first()

        if not conversation:
            raise HTTPException(status_code=404, detail="Conversation not found")

        old_status = conversation.status
        conversation.status = status

        # Assign admin if not assigned
        if not conversation.admin_id:
            conversation.admin_id = admin.firebase_uid

        db.commit()

        # Audit log
        audit = AuditService(db)
        audit.log_action(
            action="admin_chat_status_updated",
            user_id=admin.firebase_uid,
            resource_type="chat",
            resource_id=str(conversation_id),
            metadata={
                "old_status": old_status,
                "new_status": status,
                "user_id": conversation.user_id
            }
        )

        return {"status": status}
    finally:
        db.close()
