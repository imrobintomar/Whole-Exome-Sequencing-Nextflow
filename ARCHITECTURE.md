# System Architecture

## Overview

This document describes the architecture of the WES Pipeline web application.

---

## High-Level Architecture

```
┌──────────────────────────────────────────────────────────────────┐
│                         USER'S BROWSER                            │
│                                                                   │
│  ┌────────────────────────────────────────────────────────────┐  │
│  │           Frontend (React/Next.js)                         │  │
│  │                                                            │  │
│  │  • User Authentication (Login/Register)                   │  │
│  │  • File Upload Interface (Drag & Drop)                    │  │
│  │  • Job Status Dashboard                                   │  │
│  │  • Result Download Interface                              │  │
│  └────────────────────────────────────────────────────────────┘  │
│                              ▲                                    │
└──────────────────────────────┼────────────────────────────────────┘
                               │ HTTPS
                               │ (REST API + JWT Auth)
                               ▼
┌──────────────────────────────────────────────────────────────────┐
│                    VERCEL (CDN + Hosting)                         │
└──────────────────────────────────────────────────────────────────┘
                               │
                               │ HTTPS
                               │
┌──────────────────────────────┼────────────────────────────────────┐
│                  YOUR SYSTEM (Local/Server)                       │
│                              ▼                                    │
│  ┌────────────────────────────────────────────────────────────┐  │
│  │          Backend API (FastAPI + Python)                    │  │
│  │                                                            │  │
│  │  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐    │  │
│  │  │   Auth       │  │  File        │  │  Job         │    │  │
│  │  │   Manager    │  │  Manager     │  │  Manager     │    │  │
│  │  └──────┬───────┘  └──────┬───────┘  └──────┬───────┘    │  │
│  │         │                 │                 │             │  │
│  │         ▼                 ▼                 ▼             │  │
│  │  ┌─────────────────────────────────────────────────────┐ │  │
│  │  │           SQLite Database                           │ │  │
│  │  │  • Users table (auth)                               │ │  │
│  │  │  • Jobs table (tracking)                            │ │  │
│  │  └─────────────────────────────────────────────────────┘ │  │
│  └────────────────────────────────────────────────────────────┘  │
│                              │                                    │
│                              │ Subprocess                         │
│                              ▼                                    │
│  ┌────────────────────────────────────────────────────────────┐  │
│  │       Nextflow Pipeline Runner (pipeline.py)               │  │
│  │                                                            │  │
│  │  • Prepares job directory                                 │  │
│  │  • Copies FASTQ files                                     │  │
│  │  • Executes: nextflow run main.nf                         │  │
│  │  • Monitors execution                                     │  │
│  │  • Updates job status                                     │  │
│  └─────────────────────┬──────────────────────────────────────┘  │
│                        │                                          │
│                        │ Subprocess                               │
│                        ▼                                          │
│  ┌────────────────────────────────────────────────────────────┐  │
│  │         Nextflow Pipeline (main.nf)                        │  │
│  │                                                            │  │
│  │  ┌──────────┐   ┌──────────┐   ┌──────────┐              │  │
│  │  │  FastQC  │──▶│   BWA    │──▶│  GATK    │              │  │
│  │  │   (QC)   │   │  (Align) │   │  (BQSR)  │              │  │
│  │  └──────────┘   └──────────┘   └──────────┘              │  │
│  │                                     │                      │  │
│  │                                     ▼                      │  │
│  │  ┌──────────┐   ┌──────────┐   ┌──────────┐              │  │
│  │  │ Filtered │◀──│ ANNOVAR  │◀──│HaplotypeC│              │  │
│  │  │   TSV    │   │(Annotate)│   │  aller   │              │  │
│  │  └──────────┘   └──────────┘   └──────────┘              │  │
│  └────────────────────────────────────────────────────────────┘  │
│                        │                                          │
│                        │ Writes                                   │
│                        ▼                                          │
│  ┌────────────────────────────────────────────────────────────┐  │
│  │              Local Filesystem Storage                      │  │
│  │                                                            │  │
│  │  uploads/                     results/                    │  │
│  │  └── {job-id}/                └── {job-id}/               │  │
│  │      ├── R1.fastq.gz              ├── input/              │  │
│  │      └── R2.fastq.gz              └── output/             │  │
│  │                                       ├── Mapsam/         │  │
│  │                                       │   └── *.bam       │  │
│  │                                       └── Germline_VCF/   │  │
│  │                                           ├── *.vcf.gz    │  │
│  │                                           ├── *.vcf       │  │
│  │                                           └── *.tsv       │  │
│  └────────────────────────────────────────────────────────────┘  │
└──────────────────────────────────────────────────────────────────┘
```

---

## Component Details

### Frontend (Next.js + React)

**Location**: Hosted on Vercel CDN
**Technology Stack**:
- Next.js 14 (React framework)
- TypeScript
- Tailwind CSS
- Axios (HTTP client)
- React Dropzone (file upload)

**Components**:
1. **LoginForm.tsx**: User authentication
2. **RegisterForm.tsx**: New user registration
3. **Dashboard.tsx**: Main application shell
4. **UploadForm.tsx**: FASTQ file upload interface
5. **JobList.tsx**: Job monitoring and result downloads

**State Management**:
- Local state with React hooks
- JWT token stored in localStorage
- Auto-refresh jobs every 10 seconds

---

### Backend API (FastAPI + Python)

**Location**: Your local system
**Technology Stack**:
- FastAPI (async web framework)
- SQLAlchemy (ORM)
- SQLite (database)
- python-jose (JWT)
- passlib (password hashing)
- python-multipart (file uploads)

**Modules**:

#### 1. **main.py**
- FastAPI application setup
- Route definitions
- CORS middleware
- File upload handlers

#### 2. **auth.py**
- User authentication
- JWT token generation/validation
- Password hashing
- OAuth2 password flow

#### 3. **database.py**
- SQLAlchemy models (User, Job)
- Database connection
- Session management
- Enum types (JobStatus)

#### 4. **schemas.py**
- Pydantic models for request/response
- Data validation
- Serialization

#### 5. **pipeline.py**
- PipelineRunner class
- Job directory preparation
- Nextflow execution
- Status tracking
- File path resolution

#### 6. **config.py**
- Environment variable loading
- Settings management
- CORS origin parsing

---

### Database Schema

**Users Table**:
```sql
CREATE TABLE users (
    id INTEGER PRIMARY KEY,
    email VARCHAR UNIQUE NOT NULL,
    username VARCHAR UNIQUE NOT NULL,
    hashed_password VARCHAR NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
```

**Jobs Table**:
```sql
CREATE TABLE jobs (
    id INTEGER PRIMARY KEY,
    job_id VARCHAR UNIQUE NOT NULL,
    user_id INTEGER NOT NULL,
    sample_name VARCHAR NOT NULL,
    status ENUM('pending','running','completed','failed'),

    fastq_r1_path VARCHAR NOT NULL,
    fastq_r2_path VARCHAR NOT NULL,

    bam_path VARCHAR NULL,
    raw_vcf_path VARCHAR NULL,
    annotated_vcf_path VARCHAR NULL,
    filtered_tsv_path VARCHAR NULL,

    error_message VARCHAR NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    started_at TIMESTAMP NULL,
    completed_at TIMESTAMP NULL,

    FOREIGN KEY (user_id) REFERENCES users(id)
);
```

---

## Data Flow

### 1. User Registration/Login

```
Browser → POST /register → Backend → Hash password → Store in DB
Browser → POST /token → Backend → Validate credentials → Return JWT
Browser stores JWT → Include in Authorization header for all requests
```

### 2. File Upload & Job Submission

```
User selects FASTQ files
  ↓
React Dropzone validates file type
  ↓
User clicks "Submit"
  ↓
Frontend: FormData with sample_name + R1 + R2
  ↓
POST /jobs/submit (with JWT)
  ↓
Backend validates JWT → Extract user_id
  ↓
Generate unique job_id (UUID)
  ↓
Create uploads/{job_id}/ directory
  ↓
Save R1.fastq.gz and R2.fastq.gz
  ↓
Create Job record in database (status=PENDING)
  ↓
Add background task: run_job_async(job_id)
  ↓
Return job details to frontend
  ↓
Frontend switches to "My Jobs" tab
```

### 3. Pipeline Execution (Background)

```
Background task starts
  ↓
Load Job from database
  ↓
Update status → RUNNING, set started_at
  ↓
Create results/{job_id}/input/ and /output/
  ↓
Copy FASTQ files to input directory
  ↓
Build Nextflow command:
  nextflow run main.nf
    --input_dir {input}
    --output_dir {output}
    --reference {ref}
    -resume
  ↓
Execute as subprocess, capture stdout/stderr
  ↓
Wait for completion
  ↓
If success:
  - Find output files (BAM, VCF, TSV)
  - Update Job with file paths
  - Set status → COMPLETED
  ↓
If failure:
  - Capture error message
  - Set status → FAILED
  ↓
Commit to database
```

### 4. Job Monitoring

```
Frontend: GET /jobs (every 10 seconds)
  ↓
Backend validates JWT → Extract user_id
  ↓
Query: SELECT * FROM jobs WHERE user_id = ?
  ↓
Return jobs with current status
  ↓
Frontend updates UI:
  - PENDING: yellow clock icon
  - RUNNING: blue spinner icon
  - COMPLETED: green checkmark + download buttons
  - FAILED: red X + error message
```

### 5. Result Download

```
User clicks "Download BAM"
  ↓
Frontend: GET /jobs/{job_id}/download/bam
  ↓
Backend validates JWT → Check job ownership
  ↓
Load Job from database
  ↓
Check status == COMPLETED
  ↓
Retrieve bam_path from database
  ↓
Return FileResponse with file stream
  ↓
Browser downloads file
```

---

## Security Architecture

### Authentication Flow

```
1. User registers → Password hashed with bcrypt (salt + hash)
2. Stored in database: username + email + hashed_password
3. Login → Validate password → Generate JWT token
4. JWT contains: {sub: username, exp: timestamp}
5. Token signed with SECRET_KEY (HS256 algorithm)
6. Frontend stores token in localStorage
7. All API requests include: Authorization: Bearer {token}
8. Backend validates: decode JWT → check signature → check expiry
9. Extract username → load User from database
10. Verify user owns requested resources
```

### File Access Control

```
- Each job has user_id foreign key
- Download endpoint checks: job.user_id == current_user.id
- Users can only access their own uploads/results
- File paths stored in database (not user-provided)
- No directory traversal possible
```

### Network Security

```
Frontend (Vercel) ◄──HTTPS──► Backend (Your System)

- All traffic encrypted with TLS
- JWT tokens never sent in query params
- CORS restricts allowed origins
- Rate limiting recommended for production
```

---

## Scalability Considerations

### Current Limitations

1. **Single Machine**: Pipeline runs on one system
2. **Sequential Jobs**: Jobs processed one at a time
3. **Local Storage**: Files stored on local disk
4. **SQLite**: Limited concurrent writes

### Future Improvements

1. **Job Queue**: Use Celery + Redis for distributed processing
2. **Multiple Workers**: Scale horizontally with multiple systems
3. **Cloud Storage**: S3/GCS for uploads and results
4. **PostgreSQL**: Better concurrent access
5. **Load Balancer**: Distribute API requests
6. **CDN**: Serve results via CDN for faster downloads
7. **Auto-scaling**: Spin up workers based on queue depth

---

## Deployment Strategies

### Development

```
Frontend: localhost:3000 (npm run dev)
Backend: localhost:8000 (python main.py)
Database: SQLite file
Network: Local only
```

### Production (ngrok)

```
Frontend: Vercel (https://yourapp.vercel.app)
Backend: Your system exposed via ngrok (https://xxx.ngrok.io)
Database: SQLite file
Network: Internet via ngrok tunnel
```

### Production (Self-hosted)

```
Frontend: Vercel (https://yourapp.vercel.app)
Backend: Your domain (https://api.yourdomain.com)
  - Nginx reverse proxy
  - Let's Encrypt SSL
  - Port forwarding (443 → 8000)
Database: SQLite or PostgreSQL
Network: Internet via your ISP
```

### Production (Cloud)

```
Frontend: Vercel
Backend: Cloud VM (AWS EC2, DigitalOcean, etc.)
  - Docker container
  - Systemd service
  - Auto-restart
Database: PostgreSQL on RDS
Storage: S3 buckets
Network: Cloud provider's network
```

---

## Monitoring & Logging

### Application Logs

```python
# Backend logs (add logging)
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

logger.info(f"Job {job_id} started")
logger.error(f"Pipeline failed: {error}")
```

### Nextflow Logs

```
results/{job-id}/output/logs/
  ├── trace.txt       # Execution trace
  ├── timeline.html   # Timeline visualization
  ├── report.html     # Execution report
  └── pipeline.svg    # DAG visualization
```

### System Monitoring

```bash
# Disk usage
df -h
du -sh uploads/ results/

# Database size
ls -lh wes_pipeline.db

# Active jobs
sqlite3 wes_pipeline.db "SELECT * FROM jobs WHERE status='running';"

# System resources
htop
nvidia-smi  # If using GPU
```

---

## Error Handling

### Frontend Errors

```typescript
try {
  await jobApi.submitJob(...)
} catch (err) {
  if (err.response?.status === 401) {
    // Redirect to login
  } else if (err.response?.status === 413) {
    // File too large
  } else {
    // Generic error message
  }
}
```

### Backend Errors

```python
try:
  job = db.query(Job).filter(...).first()
  if not job:
    raise HTTPException(status_code=404, detail="Job not found")
except Exception as e:
  logger.error(f"Unexpected error: {e}")
  raise HTTPException(status_code=500, detail="Internal server error")
```

### Pipeline Errors

```python
# Captured in pipeline.py
if process.returncode != 0:
  job.status = JobStatus.FAILED
  job.error_message = stderr[:500]
  db.commit()
```

---

## Backup & Recovery

### Database Backup

```bash
# Daily backup script
#!/bin/bash
DATE=$(date +%Y%m%d)
cp wes_pipeline.db backups/wes_$DATE.db
find backups/ -mtime +30 -delete  # Keep 30 days
```

### File Backup

```bash
# Archive completed jobs
tar -czf archive_$DATE.tar.gz results/{job-id}/
aws s3 cp archive_$DATE.tar.gz s3://my-bucket/archives/
```

### Disaster Recovery

```bash
# Restore database
cp backups/wes_20260115.db wes_pipeline.db

# Restore files
tar -xzf archive_20260115.tar.gz -C results/
```

---

## Performance Optimization

### Database

- Add indexes on frequently queried columns
- Use connection pooling
- Vacuum SQLite periodically

### File Storage

- Use separate disk for uploads/results
- SSD for faster I/O
- Regular cleanup of old files

### API

- Enable gzip compression
- Cache static responses
- Rate limiting to prevent abuse

### Pipeline

- Use `-resume` for fault tolerance
- Optimize resource allocation
- Use SSD for work directory

---

## Maintenance Checklist

**Daily**:
- Monitor disk space
- Check for failed jobs
- Review error logs

**Weekly**:
- Clean old uploads (>30 days)
- Backup database
- Review system resources

**Monthly**:
- Update dependencies
- Security patches
- Archive old results
- Performance review

---

## Conclusion

This architecture provides a scalable, secure, and user-friendly web interface for your Whole Exome Sequencing pipeline. The modular design allows for easy upgrades and customization as your needs grow.
