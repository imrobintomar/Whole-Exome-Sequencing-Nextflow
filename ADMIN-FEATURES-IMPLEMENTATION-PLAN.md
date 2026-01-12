# Admin Dashboard Enhancement - Implementation Plan

## Phase 1: Core Improvements (Week 1-2)

### 1.1 Real-time Monitoring
**Files to modify:**
- `frontend/components/AdminDashboard.tsx`

**Implementation:**
```typescript
// Add state for auto-refresh
const [autoRefresh, setAutoRefresh] = useState(true);
const [refreshInterval, setRefreshInterval] = useState(30000); // 30 seconds

useEffect(() => {
  if (!autoRefresh) return;

  const interval = setInterval(() => {
    fetchDashboardStats();
    fetchJobs();
    fetchSystemMetrics();
  }, refreshInterval);

  return () => clearInterval(interval);
}, [autoRefresh, refreshInterval]);

// Add UI controls
<div className="flex items-center gap-4">
  <button onClick={() => setAutoRefresh(!autoRefresh)}>
    {autoRefresh ? '⏸️ Pause' : '▶️ Resume'} Auto-refresh
  </button>
  <select value={refreshInterval} onChange={(e) => setRefreshInterval(Number(e.target.value))}>
    <option value={10000}>10s</option>
    <option value={30000}>30s</option>
    <option value={60000}>1m</option>
  </select>
</div>
```

### 1.2 Job Filtering & Search
**Backend changes:**
- `backend/modules/admin/routes.py` - Already has basic filtering

**Frontend changes:**
- `frontend/components/AdminDashboard.tsx`

**Implementation:**
```typescript
// Add filter state
const [filters, setFilters] = useState({
  status: '',
  userId: '',
  dateFrom: '',
  dateTo: '',
  searchQuery: ''
});

// Add search/filter UI
<div className="filters-panel">
  <input
    type="text"
    placeholder="Search by job ID or sample name..."
    value={filters.searchQuery}
    onChange={(e) => setFilters({...filters, searchQuery: e.target.value})}
  />

  <select value={filters.status} onChange={(e) => setFilters({...filters, status: e.target.value})}>
    <option value="">All Status</option>
    <option value="pending">Pending</option>
    <option value="running">Running</option>
    <option value="completed">Completed</option>
    <option value="failed">Failed</option>
  </select>

  <input
    type="date"
    placeholder="From"
    value={filters.dateFrom}
    onChange={(e) => setFilters({...filters, dateFrom: e.target.value})}
  />

  <input
    type="date"
    placeholder="To"
    value={filters.dateTo}
    onChange={(e) => setFilters({...filters, dateTo: e.target.value})}
  />
</div>

// Apply filters
const filteredJobs = jobs.filter(job => {
  if (filters.status && job.status !== filters.status) return false;
  if (filters.searchQuery && !job.job_id.includes(filters.searchQuery) && !job.sample_name.includes(filters.searchQuery)) return false;
  if (filters.dateFrom && new Date(job.created_at) < new Date(filters.dateFrom)) return false;
  if (filters.dateTo && new Date(job.created_at) > new Date(filters.dateTo)) return false;
  return true;
});
```

### 1.3 Bulk Job Actions
**Backend additions:**
- `backend/modules/admin/routes.py`

**New endpoint:**
```python
@router.post("/jobs/bulk-action")
async def bulk_job_action(
    action: str,  # 'cancel', 'delete', 'rerun'
    job_ids: list[str],
    admin=Depends(require_admin)
):
    """
    Perform bulk actions on multiple jobs
    """
    db = SessionLocal()
    try:
        results = {
            "success": [],
            "failed": []
        }

        for job_id in job_ids:
            try:
                job = db.query(Job).filter(Job.job_id == job_id).first()
                if not job:
                    results["failed"].append({"job_id": job_id, "error": "Not found"})
                    continue

                if action == "delete":
                    db.delete(job)
                    results["success"].append(job_id)
                elif action == "cancel":
                    # Kill Nextflow process
                    PipelineRunner.cancel_job(job_id)
                    job.status = "cancelled"
                    results["success"].append(job_id)
                elif action == "rerun":
                    # Implement rerun logic
                    results["success"].append(job_id)

            except Exception as e:
                results["failed"].append({"job_id": job_id, "error": str(e)})

        db.commit()
        return results
    finally:
        db.close()
```

**Frontend:**
```typescript
// Add selection state
const [selectedJobs, setSelectedJobs] = useState<string[]>([]);

// Checkbox for each job
<input
  type="checkbox"
  checked={selectedJobs.includes(job.job_id)}
  onChange={(e) => {
    if (e.target.checked) {
      setSelectedJobs([...selectedJobs, job.job_id]);
    } else {
      setSelectedJobs(selectedJobs.filter(id => id !== job.job_id));
    }
  }}
/>

// Bulk action buttons
<div className="bulk-actions">
  <button onClick={() => bulkAction('cancel')}>Cancel Selected ({selectedJobs.length})</button>
  <button onClick={() => bulkAction('delete')}>Delete Selected ({selectedJobs.length})</button>
</div>

const bulkAction = async (action: string) => {
  if (!confirm(`${action} ${selectedJobs.length} jobs?`)) return;

  const response = await api.post('/admin/jobs/bulk-action', {
    action,
    job_ids: selectedJobs
  });

  alert(`Success: ${response.success.length}, Failed: ${response.failed.length}`);
  setSelectedJobs([]);
  fetchJobs();
};
```

---

## Phase 2: User Management (Week 3)

### 2.1 User Actions (Ban/Suspend)
**Backend:**
```python
# Add to User model in database.py
class User(Base):
    # ... existing fields
    is_active = Column(Boolean, default=True)
    is_banned = Column(Boolean, default=False)
    ban_reason = Column(String, nullable=True)
    banned_at = Column(DateTime, nullable=True)

# New endpoint in backend/modules/admin/routes.py
@router.post("/users/{user_uid}/ban")
async def ban_user(
    user_uid: str,
    reason: str,
    admin=Depends(require_admin)
):
    db = SessionLocal()
    try:
        user = db.query(User).filter(User.firebase_uid == user_uid).first()
        if not user:
            raise HTTPException(status_code=404, detail="User not found")

        user.is_banned = True
        user.ban_reason = reason
        user.banned_at = datetime.now()
        db.commit()

        # Log action
        AuditService.log_action(admin.firebase_uid, "ban_user", user_uid, reason)

        return {"success": True, "message": f"User {user.email} banned"}
    finally:
        db.close()

@router.post("/users/{user_uid}/unban")
async def unban_user(user_uid: str, admin=Depends(require_admin)):
    # Similar implementation
    pass
```

**Middleware update:**
```python
# Update firebase_auth.py to check ban status
async def get_current_user(...):
    # ... existing code

    if user.is_banned:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail=f"Account banned: {user.ban_reason}",
        )

    return user
```

### 2.2 Export User List
**Backend:**
```python
@router.get("/users/export")
async def export_users(admin=Depends(require_admin)):
    """
    Export all users as CSV
    """
    db = SessionLocal()
    try:
        users = db.query(User).all()

        csv_data = "UID,Email,Created At,Status,Jobs Count\n"
        for user in users:
            jobs_count = db.query(Job).filter(Job.user_id == user.firebase_uid).count()
            csv_data += f"{user.firebase_uid},{user.email},{user.created_at},{user.is_banned and 'Banned' or 'Active'},{jobs_count}\n"

        return Response(
            content=csv_data,
            media_type="text/csv",
            headers={"Content-Disposition": "attachment; filename=users.csv"}
        )
    finally:
        db.close()
```

**Frontend:**
```typescript
const exportUsers = async () => {
  const response = await api.get('/admin/users/export', {
    responseType: 'blob'
  });

  const url = window.URL.createObjectURL(new Blob([response.data]));
  const link = document.createElement('a');
  link.href = url;
  link.download = `users_${new Date().toISOString().split('T')[0]}.csv`;
  link.click();
};
```

---

## Phase 3: System Health & Alerts (Week 4)

### 3.1 Health Monitoring Service
**New file: `backend/services/health_monitor.py`**
```python
import psutil
from datetime import datetime
from typing import Dict, List
import smtplib
from email.mime.text import MIMEText

class HealthMonitor:
    THRESHOLDS = {
        "cpu": 90,
        "memory": 90,
        "disk": 90
    }

    @staticmethod
    def check_system_health() -> Dict:
        """
        Check if system metrics exceed thresholds
        """
        alerts = []

        cpu = psutil.cpu_percent(interval=1)
        memory = psutil.virtual_memory().percent
        disk = psutil.disk_usage('/').percent

        if cpu > HealthMonitor.THRESHOLDS["cpu"]:
            alerts.append({
                "type": "cpu",
                "value": cpu,
                "threshold": HealthMonitor.THRESHOLDS["cpu"],
                "severity": "critical",
                "message": f"CPU usage at {cpu}%"
            })

        if memory > HealthMonitor.THRESHOLDS["memory"]:
            alerts.append({
                "type": "memory",
                "value": memory,
                "threshold": HealthMonitor.THRESHOLDS["memory"],
                "severity": "critical",
                "message": f"Memory usage at {memory}%"
            })

        if disk > HealthMonitor.THRESHOLDS["disk"]:
            alerts.append({
                "type": "disk",
                "value": disk,
                "threshold": HealthMonitor.THRESHOLDS["disk"],
                "severity": "warning",
                "message": f"Disk usage at {disk}%"
            })

        return {
            "status": "healthy" if len(alerts) == 0 else "unhealthy",
            "alerts": alerts,
            "checked_at": datetime.now().isoformat()
        }

    @staticmethod
    def send_alert_email(alert: Dict):
        """
        Send email notification for critical alerts
        """
        # Implement email sending
        pass
```

**Backend endpoint:**
```python
@router.get("/health/check")
async def health_check(admin=Depends(require_admin)):
    """
    Get system health status
    """
    return HealthMonitor.check_system_health()
```

**Frontend:**
```typescript
// Add to AdminDashboard
const [healthAlerts, setHealthAlerts] = useState([]);

useEffect(() => {
  const checkHealth = async () => {
    const health = await api.get('/admin/health/check');
    if (health.data.alerts.length > 0) {
      setHealthAlerts(health.data.alerts);
      // Show notification
      showNotification('System Health Alert', health.data.alerts[0].message);
    }
  };

  // Check every 5 minutes
  const interval = setInterval(checkHealth, 5 * 60 * 1000);
  return () => clearInterval(interval);
}, []);

// Display alerts
{healthAlerts.length > 0 && (
  <div className="health-alerts">
    {healthAlerts.map(alert => (
      <div className={`alert alert-${alert.severity}`}>
        ⚠️ {alert.message}
      </div>
    ))}
  </div>
)}
```

---

## Phase 4: Analytics & Reports (Week 5-6)

### 4.1 Usage Trends Backend
**New file: `backend/services/analytics_service.py`**
```python
from sqlalchemy import func, extract
from datetime import datetime, timedelta

class AnalyticsService:
    @staticmethod
    def get_jobs_trend(days: int = 30):
        """
        Get job creation trend for last N days
        """
        db = SessionLocal()
        try:
            start_date = datetime.now() - timedelta(days=days)

            results = db.query(
                func.date(Job.created_at).label('date'),
                func.count(Job.id).label('count')
            ).filter(
                Job.created_at >= start_date
            ).group_by(
                func.date(Job.created_at)
            ).order_by('date').all()

            return [
                {"date": str(r.date), "count": r.count}
                for r in results
            ]
        finally:
            db.close()

    @staticmethod
    def get_revenue_trend(months: int = 12):
        """
        Monthly revenue for last N months
        """
        db = SessionLocal()
        try:
            start_date = datetime.now() - timedelta(days=months * 30)

            results = db.query(
                extract('year', Subscription.created_at).label('year'),
                extract('month', Subscription.created_at).label('month'),
                func.sum(SubscriptionPlan.price_cents).label('revenue')
            ).join(
                SubscriptionPlan
            ).filter(
                Subscription.status == 'active',
                Subscription.created_at >= start_date
            ).group_by('year', 'month').order_by('year', 'month').all()

            return [
                {
                    "month": f"{int(r.year)}-{int(r.month):02d}",
                    "revenue_usd": r.revenue / 100
                }
                for r in results
            ]
        finally:
            db.close()

    @staticmethod
    def get_user_retention():
        """
        Calculate user retention rate
        """
        db = SessionLocal()
        try:
            total_users = db.query(User).count()

            # Active in last 30 days
            thirty_days_ago = datetime.now() - timedelta(days=30)
            active_users = db.query(User).join(Job).filter(
                Job.created_at >= thirty_days_ago
            ).distinct().count()

            retention_rate = (active_users / total_users * 100) if total_users > 0 else 0

            return {
                "total_users": total_users,
                "active_users": active_users,
                "retention_rate": round(retention_rate, 2)
            }
        finally:
            db.close()
```

**Backend endpoints:**
```python
@router.get("/analytics/jobs-trend")
async def get_jobs_trend(days: int = 30, admin=Depends(require_admin)):
    return AnalyticsService.get_jobs_trend(days)

@router.get("/analytics/revenue-trend")
async def get_revenue_trend(months: int = 12, admin=Depends(require_admin)):
    return AnalyticsService.get_revenue_trend(months)

@router.get("/analytics/retention")
async def get_retention(admin=Depends(require_admin)):
    return AnalyticsService.get_user_retention()
```

### 4.2 Charts Frontend
**Install chart library:**
```bash
npm install recharts
```

**Frontend component:**
```typescript
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend } from 'recharts';

const JobsTrendChart = () => {
  const [data, setData] = useState([]);

  useEffect(() => {
    const fetchTrend = async () => {
      const response = await api.get('/admin/analytics/jobs-trend?days=30');
      setData(response.data);
    };
    fetchTrend();
  }, []);

  return (
    <div className="chart-container">
      <h3>Jobs Created (Last 30 Days)</h3>
      <LineChart width={600} height={300} data={data}>
        <CartesianGrid strokeDasharray="3 3" />
        <XAxis dataKey="date" />
        <YAxis />
        <Tooltip />
        <Legend />
        <Line type="monotone" dataKey="count" stroke="#8884d8" />
      </LineChart>
    </div>
  );
};
```

---

## Phase 5: Advanced Features (Week 7-8)

### 5.1 Activity Log
**New table in `database_extensions.py`:**
```python
class ActivityLog(Base):
    __tablename__ = "activity_logs"

    id = Column(Integer, primary_key=True, index=True)
    user_id = Column(String, index=True)
    admin_id = Column(String, nullable=True)  # If admin action
    action = Column(String)  # 'login', 'job_submit', 'subscription_change', etc.
    details = Column(JSON)
    ip_address = Column(String, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow)
```

**Service:**
```python
class ActivityLogService:
    @staticmethod
    def log(user_id: str, action: str, details: dict = None, ip: str = None):
        db = SessionLocal()
        try:
            log_entry = ActivityLog(
                user_id=user_id,
                action=action,
                details=details,
                ip_address=ip
            )
            db.add(log_entry)
            db.commit()
        finally:
            db.close()
```

### 5.2 Notification System
**New file: `backend/services/notification_service.py`:**
```python
class NotificationService:
    @staticmethod
    def send_email(to: str, subject: str, body: str):
        # Implement SMTP
        pass

    @staticmethod
    def send_slack_notification(message: str, webhook_url: str):
        import requests
        requests.post(webhook_url, json={"text": message})

    @staticmethod
    def notify_quota_warning(user_email: str, usage_percent: float):
        subject = "Usage Quota Warning"
        body = f"You've used {usage_percent}% of your monthly quota."
        NotificationService.send_email(user_email, subject, body)
```

---

## Implementation Priority

### Week 1-2: Quick Wins
1. ✅ Real-time auto-refresh
2. ✅ Job filtering & search
3. ✅ Bulk job actions
4. ✅ Export users CSV

### Week 3-4: Core Features
5. ✅ User ban/suspend
6. ✅ System health monitoring
7. ✅ Health alerts

### Week 5-6: Analytics
8. ✅ Usage trends charts
9. ✅ Revenue analytics
10. ✅ User retention metrics

### Week 7-8: Advanced
11. ✅ Activity logging
12. ✅ Notification system
13. ✅ Data cleanup tools

---

## Files to Create/Modify

### Backend
- ✅ `backend/services/health_monitor.py` (NEW)
- ✅ `backend/services/analytics_service.py` (NEW)
- ✅ `backend/services/notification_service.py` (NEW)
- ✅ `backend/modules/admin/routes.py` (MODIFY - add endpoints)
- ✅ `backend/database_extensions.py` (MODIFY - add ActivityLog)

### Frontend
- ✅ `frontend/components/AdminDashboard.tsx` (MODIFY - main updates)
- ✅ `frontend/components/JobsTrendChart.tsx` (NEW)
- ✅ `frontend/components/RevenueDashboard.tsx` (NEW)
- ✅ `frontend/components/HealthAlerts.tsx` (NEW)
- ✅ `frontend/lib/api.ts` (MODIFY - add new endpoints)

### Database Migrations
- ✅ Add `is_banned`, `ban_reason`, `banned_at` to User
- ✅ Create `activity_logs` table

---

## Testing Checklist

- [ ] Auto-refresh works and can be paused
- [ ] Job filtering by all criteria
- [ ] Bulk actions on multiple jobs
- [ ] User export downloads CSV
- [ ] User ban prevents login
- [ ] Health alerts trigger at thresholds
- [ ] Charts render correctly with data
- [ ] Activity logs capture all actions
- [ ] Email notifications send

---

## Environment Variables Needed

```env
# Email settings
SMTP_HOST=smtp.gmail.com
SMTP_PORT=587
SMTP_USER=your-email@gmail.com
SMTP_PASSWORD=your-app-password
ADMIN_EMAIL=admin@atgcflow.com

# Slack (optional)
SLACK_WEBHOOK_URL=https://hooks.slack.com/services/YOUR/WEBHOOK/URL

# Health thresholds
HEALTH_CPU_THRESHOLD=90
HEALTH_MEMORY_THRESHOLD=90
HEALTH_DISK_THRESHOLD=90
```

---

This plan provides a structured approach to implementing all requested features systematically.
