# SaaS Extension - Implementation Documentation

## Overview

This document describes the SaaS features that have been **additively** added to the WES Platform backend without modifying any existing code.

## Architecture

### Design Principles
- **Additive Only**: Zero modifications to existing code
- **Environment-Driven**: All sensitive config via ENV variables
- **Graceful Degradation**: Platform works with or without SaaS modules
- **Security First**: Admin guards, audit logging, webhook verification

### Module Structure

```
backend/
├── database_extensions.py          # New database tables (8 tables)
├── init_saas_db.py                 # Database initialization script
├── .env.saas.example               # Environment template
│
├── middleware/
│   ├── admin_guard.py              # Admin-only access control
│   └── subscription_guard.py       # Usage limit enforcement
│
├── services/
│   ├── audit_service.py            # Audit logging
│   ├── usage_tracking_service.py   # Monthly usage tracking
│   ├── stripe_service.py           # Stripe API wrapper
│   └── metrics_service.py          # System metrics (CPU/memory/disk)
│
└── modules/
    ├── admin/
    │   └── routes.py               # Admin API endpoints
    ├── billing/
    │   └── routes.py               # Billing & subscription endpoints
    └── chat/
        └── routes.py               # User-Admin chat system
```

## Database Schema

### New Tables (8 total)

1. **subscription_plans** - Plan definitions (Free, Basic, Pro)
2. **subscriptions** - User subscriptions (linked to Stripe)
3. **usage_tracking** - Monthly job usage counters
4. **chat_conversations** - Support ticket threads
5. **chat_messages** - Messages within conversations
6. **admin_users** - Admin user registry
7. **audit_log** - Compliance and security audit trail
8. **webhook_events** - Stripe webhook idempotency

## API Endpoints

### Admin Routes (`/admin/*`)
All require admin authentication via `ADMIN_USER_UIDS` env variable.

- `GET /admin/dashboard/stats` - Dashboard statistics
- `GET /admin/jobs` - List all jobs (filterable)
- `GET /admin/jobs/{job_id}/logs` - View Nextflow logs
- `GET /admin/users` - List all users with subscription info
- `GET /admin/system/metrics` - Real-time system metrics

### Billing Routes (`/billing/*`)

**Public:**
- `GET /billing/plans` - View subscription plans

**User (authenticated):**
- `GET /billing/subscription` - Current subscription details
- `POST /billing/checkout` - Create Stripe checkout session
- `POST /billing/portal` - Access customer portal
- `GET /billing/usage` - Current month usage stats

**Webhook (public, verified):**
- `POST /billing/webhooks/stripe` - Stripe event handler

### Chat Routes (`/chat/*`)

**User:**
- `POST /chat/conversations` - Create new conversation
- `GET /chat/conversations` - List my conversations
- `GET /chat/conversations/{id}/messages` - View messages
- `POST /chat/conversations/{id}/messages` - Send message
- `PATCH /chat/conversations/{id}/close` - Close conversation

**Admin:**
- `GET /chat/admin/conversations` - List all conversations
- `GET /chat/admin/conversations/{id}/messages` - View any conversation
- `POST /chat/admin/conversations/{id}/messages` - Reply as admin
- `PATCH /chat/admin/conversations/{id}/status` - Update status

### Job Submission Wrapper

**New Endpoint:**
- `POST /jobs/submit-with-billing` - Submit job with usage enforcement

**Original Endpoint (unchanged):**
- `POST /jobs/submit` - Original submission (backward compatible)

## Setup Instructions

### 1. Install Dependencies

```bash
cd backend
pip install -r requirements.txt
```

New dependencies:
- `stripe>=10.0.0` - Stripe payment processing
- `psutil>=5.9.0` - System metrics

### 2. Configure Environment

```bash
cp .env.saas.example .env
```

**Required variables:**
```bash
FRONTEND_URL=http://localhost:3000
STRIPE_SECRET_KEY=sk_test_xxxxx
STRIPE_PUBLISHABLE_KEY=pk_test_xxxxx
STRIPE_WEBHOOK_SECRET=whsec_xxxxx
ADMIN_USER_UIDS=firebase_uid_1,firebase_uid_2
```

### 3. Initialize Database

```bash
python3 init_saas_db.py
```

This creates:
- All 8 new tables
- 3 default subscription plans (Free, Basic, Pro)

### 4. Configure Stripe

#### A. Create Products in Stripe Dashboard

1. Go to https://dashboard.stripe.com/products
2. Create products:
   - **Basic**: $29/month, 10 jobs/month
   - **Pro**: $99/month, 50 jobs/month
3. Copy the Price IDs (e.g., `price_1234567890abcdef`)

#### B. Update Database

```sql
UPDATE subscription_plans
SET stripe_price_id = 'price_xxxxx'
WHERE name = 'Basic';

UPDATE subscription_plans
SET stripe_price_id = 'price_yyyyy'
WHERE name = 'Pro';
```

#### C. Create Webhook

1. Go to https://dashboard.stripe.com/webhooks
2. Add endpoint: `https://your-backend.com/billing/webhooks/stripe`
3. Select events:
   - `checkout.session.completed`
   - `customer.subscription.updated`
   - `customer.subscription.deleted`
   - `invoice.payment_failed`
4. Copy webhook secret to `STRIPE_WEBHOOK_SECRET`

### 5. Add Admin Users

**Option A: Via Database**
```sql
INSERT INTO admin_users (firebase_uid, email)
VALUES ('your_firebase_uid', 'admin@example.com');
```

**Option B: Via ENV Variable**
```bash
ADMIN_USER_UIDS=uid1,uid2,uid3
```

### 6. Start Backend

```bash
uvicorn main:app --reload
```

Look for: `✅ SaaS modules loaded successfully`

## Testing

### 1. Verify Module Loading

```bash
python3 -c "from modules.admin.routes import router; print('✓ Admin routes')"
python3 -c "from modules.billing.routes import router; print('✓ Billing routes')"
python3 -c "from modules.chat.routes import router; print('✓ Chat routes')"
```

### 2. Test Endpoints

```bash
# View subscription plans (public)
curl http://localhost:8000/billing/plans

# Admin dashboard (requires admin auth)
curl -H "Authorization: Bearer YOUR_FIREBASE_TOKEN" \
     http://localhost:8000/admin/dashboard/stats

# Current subscription (requires user auth)
curl -H "Authorization: Bearer YOUR_FIREBASE_TOKEN" \
     http://localhost:8000/billing/subscription
```

### 3. Test Stripe Integration

Use Stripe test cards:
- Success: `4242 4242 4242 4242`
- Decline: `4000 0000 0000 0002`

1. Call `/billing/checkout` with plan_id
2. Complete checkout with test card
3. Verify webhook received at `/billing/webhooks/stripe`
4. Check subscription created in database

### 4. Test Usage Limits

```bash
# Submit job with billing enforcement
curl -X POST http://localhost:8000/jobs/submit-with-billing \
  -H "Authorization: Bearer YOUR_TOKEN" \
  -F "sample_name=test" \
  -F "fastq_r1=@test_R1.fastq.gz" \
  -F "fastq_r2=@test_R2.fastq.gz"
```

Expected behavior:
- Free tier: 2 jobs/month
- After limit: HTTP 429 "Limit reached"
- With subscription: Higher limit

## Features

### ✅ Completed (Phase 1 - MVP)

1. **Subscription Management**
   - Free tier (2 jobs/month)
   - Paid plans (Basic: 10 jobs, Pro: 50 jobs)
   - Stripe checkout integration
   - Customer portal for self-service
   - Automatic usage tracking

2. **Admin Dashboard**
   - User statistics
   - Job monitoring
   - Revenue metrics
   - System health (CPU/memory/disk)
   - Nextflow process count
   - User management

3. **Chat Support**
   - User can create conversations
   - Requires active subscription with chat support
   - Admin can view/reply to all conversations
   - Status tracking (open/resolved/closed)

4. **Usage Enforcement**
   - Monthly job limits per plan
   - Row-level locking for concurrency
   - Graceful error messages
   - Automatic counter reset each month

5. **Audit Logging**
   - All sensitive operations logged
   - User ID, action, resource, metadata
   - Compliance and security tracking

6. **Webhook Handling**
   - Stripe event processing
   - Idempotency (no duplicate processing)
   - Automatic subscription updates
   - Payment failure handling

### 🔄 Pending (Phase 2-4)

**Phase 2: Enhanced Chat (Week 3-4)**
- Real-time updates (WebSockets/Server-Sent Events)
- Email notifications for new messages
- Chat analytics dashboard
- Canned responses for admins

**Phase 3: Monitoring & Analytics (Week 5-6)**
- Usage analytics dashboard
- Performance metrics visualization
- User behavior tracking
- Revenue forecasting

**Phase 4: Production Hardening (Week 7-8)**
- Rate limiting (Redis/memory-based)
- Database backups (automated)
- Security audit (OWASP compliance)
- Load testing (locust/k6)
- Documentation (API docs, user guides)

## Security

### Implemented

✅ **Authentication**
- Firebase JWT verification
- Admin guard middleware
- User-resource ownership checks

✅ **Authorization**
- Admin-only routes protected
- User can only access own resources
- Subscription requirements enforced

✅ **Data Protection**
- All secrets via environment variables
- Stripe webhook signature verification
- SQL injection protection (SQLAlchemy ORM)

✅ **Audit Trail**
- All sensitive actions logged
- User ID, timestamp, metadata captured
- Immutable audit log

### Production Recommendations

🔒 **HTTPS Required**
- Use Let's Encrypt or similar
- Enforce HTTPS redirects

🔒 **Rate Limiting**
- Implement per-user rate limits
- Protect webhook endpoints

🔒 **Database Security**
- Switch to PostgreSQL for production
- Enable SSL/TLS connections
- Regular backups

🔒 **Monitoring**
- Set up error tracking (Sentry)
- Log aggregation (ELK/Datadog)
- Uptime monitoring

## Troubleshooting

### SaaS Modules Not Loading

**Error:** `⚠️ SaaS modules not available`

**Solution:**
1. Check file structure: `ls modules/admin/routes.py`
2. Install dependencies: `pip install -r requirements.txt`
3. Check for import errors: `python3 -c "from modules.admin.routes import router"`

### Stripe Webhooks Failing

**Error:** `Invalid signature`

**Solution:**
1. Verify `STRIPE_WEBHOOK_SECRET` matches Stripe dashboard
2. Check webhook endpoint URL is correct
3. Ensure raw body is sent (FastAPI does this automatically)

### Usage Limits Not Enforced

**Error:** Users can submit unlimited jobs

**Solution:**
1. Frontend must call `/jobs/submit-with-billing` (not `/jobs/submit`)
2. Verify `subscription_guard` middleware is working
3. Check `usage_tracking` table for correct counts

### Admin Dashboard Returns 403

**Error:** `Not authorized`

**Solution:**
1. Add Firebase UID to `ADMIN_USER_UIDS` in .env
2. OR insert into `admin_users` table
3. Verify Firebase token is valid

## Migration Path

If you have existing users, follow this migration:

### 1. Preserve Existing Data
```sql
-- Backup existing database
sqlite3 wes_pipeline.db ".backup backup.db"
```

### 2. Add SaaS Tables
```bash
python3 init_saas_db.py
```

### 3. Assign Default Plans
```sql
-- Create usage records for existing users with free tier
INSERT INTO usage_tracking (user_id, month, jobs_executed, jobs_limit)
SELECT uid, strftime('%Y%m', 'now'), 0, 2
FROM users;
```

### 4. No Breaking Changes
- Existing `/jobs/submit` endpoint still works
- Users continue normal operations
- Admins can gradually migrate users to `/jobs/submit-with-billing`

## Environment Variables Reference

| Variable | Required | Default | Description |
|----------|----------|---------|-------------|
| `FRONTEND_URL` | Yes | - | Frontend URL for Stripe redirects |
| `STRIPE_SECRET_KEY` | Yes | - | Stripe API secret key |
| `STRIPE_PUBLISHABLE_KEY` | No | - | Stripe publishable key (for frontend) |
| `STRIPE_WEBHOOK_SECRET` | Yes | - | Stripe webhook signing secret |
| `ADMIN_USER_UIDS` | Yes | - | Comma-separated Firebase UIDs for admins |
| `DATABASE_URL` | No | sqlite | Database connection string |
| `DEFAULT_FREE_JOBS_LIMIT` | No | 2 | Free tier monthly job limit |
| `SESSION_TIMEOUT_HOURS` | No | 24 | Session timeout duration |

## Support

For issues or questions:
1. Check this documentation first
2. Review `.env.saas.example` for configuration
3. Check audit logs: `SELECT * FROM audit_log ORDER BY created_at DESC LIMIT 20`
4. Review Stripe webhook logs in dashboard

## License

This SaaS extension follows the same license as the main WES Platform project.
