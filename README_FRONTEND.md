# WES Pipeline with Web Frontend

Complete Whole Exome Sequencing pipeline with web-based frontend for remote job submission and result download.

## Architecture

```
┌─────────────────┐         ┌──────────────────┐         ┌──────────────────┐
│   Frontend      │         │   Backend API    │         │    Nextflow      │
│   (Vercel)      │ ◄────── │   (Your System)  │ ◄────── │    Pipeline      │
│   React/Next.js │  HTTPS  │   FastAPI        │  Local  │   (Your System)  │
└─────────────────┘         └──────────────────┘         └──────────────────┘
        │                            │                             │
        │                            │                             │
        └────── User uploads ────────┤                             │
                FASTQ files          │                             │
                                     ├──── Stores files ───────────┤
                                     │     & runs pipeline         │
                                     │                             │
                                     ├──── Tracks job status ──────┤
                                     │                             │
        ┌──── Downloads results ─────┤                             │
        │      (BAM, VCF, TSV)       │                             │
        │                            │                             │
```

## Features

- User authentication (registration & login)
- Paired-end FASTQ file upload via web interface
- Automatic pipeline execution on your system
- Real-time job status tracking
- Download results: BAM, Raw VCF, Annotated VCF, Filtered TSV
- SQLite database for job tracking
- Responsive UI with drag-and-drop file upload

---

## Setup Instructions

### 1. Backend Setup (Your System)

The backend API runs on your local system and manages file uploads, pipeline execution, and results.

#### Step 1: Install Python Dependencies

```bash
cd backend
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -r requirements.txt
```

#### Step 2: Configure Environment

```bash
cp .env.example .env
nano .env  # Edit with your configuration
```

Update the `.env` file with your paths:

```bash
SECRET_KEY=generate-a-random-secret-key-here
ALGORITHM=HS256
ACCESS_TOKEN_EXPIRE_MINUTES=30

# Storage paths
UPLOAD_DIR=/media/drprabudh/m3/Nextflow-Script/WholeExome/uploads
RESULTS_DIR=/media/drprabudh/m3/Nextflow-Script/WholeExome/results

# Nextflow configuration
NEXTFLOW_SCRIPT=/media/drprabudh/m3/Nextflow-Script/WholeExome/main.nf
REFERENCE_GENOME=/media/drprabudh/m1/hg38/hg38.fa
KNOWN_SITES_1=/media/drprabudh/m1/vcf_file/Homo_sapiens_assembly38.known_indels.vcf.gz
KNOWN_SITES_2=/media/drprabudh/m1/vcf_file/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

# Database
DATABASE_URL=sqlite:///./wes_pipeline.db

# CORS - Add your Vercel URL here after deployment
CORS_ORIGINS=http://localhost:3000,https://your-app.vercel.app
```

#### Step 3: Run Backend Server

```bash
cd backend
source venv/bin/activate
python main.py
```

The API will be available at `http://localhost:8000`

#### Step 4: Expose Backend to Internet

Since your backend is on your local system and frontend is on Vercel, you need to expose your backend. Options:

**Option A: ngrok (Easiest for testing)**
```bash
# Install ngrok: https://ngrok.com/download
ngrok http 8000
```
Copy the HTTPS URL (e.g., `https://abc123.ngrok.io`)

**Option B: Your own domain with port forwarding**
- Configure your router to forward port 8000
- Use a dynamic DNS service (e.g., No-IP, DuckDNS)
- Set up SSL with Let's Encrypt

**Option C: Tailscale (Secure)**
```bash
# Install Tailscale on your system
# Frontend can connect via Tailscale network
```

---

### 2. Frontend Setup (Vercel)

#### Step 1: Install Dependencies Locally (for testing)

```bash
cd frontend
npm install
```

#### Step 2: Configure Environment

```bash
cp .env.local.example .env.local
nano .env.local
```

Update with your backend URL:
```bash
NEXT_PUBLIC_API_URL=http://localhost:8000  # For local testing
# OR
NEXT_PUBLIC_API_URL=https://your-ngrok-url.ngrok.io  # For production
```

#### Step 3: Test Locally

```bash
npm run dev
```

Visit `http://localhost:3000` to test the frontend.

#### Step 4: Deploy to Vercel

1. Push your code to GitHub:
```bash
git init
git add frontend/
git commit -m "Add WES pipeline frontend"
git push origin main
```

2. Go to [vercel.com](https://vercel.com) and sign in

3. Click "New Project" and import your GitHub repository

4. Configure project:
   - **Framework Preset**: Next.js
   - **Root Directory**: `frontend`
   - **Build Command**: `npm run build`
   - **Output Directory**: `.next`

5. Add Environment Variable:
   - Key: `NEXT_PUBLIC_API_URL`
   - Value: `https://your-backend-url` (ngrok or your domain)

6. Click "Deploy"

7. Update Backend CORS:
   - After deployment, copy your Vercel URL (e.g., `https://your-app.vercel.app`)
   - Add it to `backend/.env` in `CORS_ORIGINS`
   - Restart backend server

---

## Usage

### For Users

1. **Register an Account**
   - Visit your Vercel URL
   - Click "Register here"
   - Enter email, username, and password

2. **Login**
   - Enter your credentials

3. **Submit Analysis**
   - Go to "Upload Files" tab
   - Enter sample name
   - Drag & drop or select R1 FASTQ file
   - Drag & drop or select R2 FASTQ file
   - Click "Submit Analysis"

4. **Monitor Progress**
   - Switch to "My Jobs" tab
   - Jobs show status: PENDING → RUNNING → COMPLETED
   - Page auto-refreshes every 10 seconds

5. **Download Results**
   - When job is COMPLETED, download buttons appear
   - Download: BAM, Raw VCF, Annotated VCF, Filtered TSV

---

## System Requirements

### Backend System
- **OS**: Linux (Ubuntu/Debian recommended)
- **CPU**: 16+ cores recommended
- **RAM**: 64GB+ recommended
- **Storage**: 500GB+ for results
- **Software**:
  - Python 3.9+
  - Nextflow
  - Docker (for pipeline containers)
  - All pipeline dependencies (BWA, GATK, ANNOVAR, etc.)

### Network
- Stable internet connection
- Ability to expose port 8000 (or use ngrok)
- Sufficient bandwidth for FASTQ uploads

---

## API Documentation

Once backend is running, visit:
- Interactive API docs: `http://localhost:8000/docs`
- Alternative docs: `http://localhost:8000/redoc`

### Key Endpoints

- `POST /register` - Create user account
- `POST /token` - Login and get JWT token
- `GET /users/me` - Get current user info
- `POST /jobs/submit` - Submit new analysis job
- `GET /jobs` - List all jobs for current user
- `GET /jobs/{job_id}` - Get job details
- `GET /jobs/{job_id}/download/{file_type}` - Download result file

---

## Security Considerations

1. **Change SECRET_KEY**: Generate a secure random key for production
   ```bash
   python -c "import secrets; print(secrets.token_urlsafe(32))"
   ```

2. **Use HTTPS**: Always use HTTPS in production (ngrok provides this)

3. **Firewall**: Only expose port 8000, block all other ports

4. **File Validation**: Backend validates FASTQ file extensions

5. **User Isolation**: Each user can only access their own jobs

6. **Token Expiration**: JWT tokens expire after 30 minutes

---

## Troubleshooting

### Backend won't start
- Check Python version: `python3 --version` (need 3.9+)
- Verify all dependencies installed: `pip list`
- Check port 8000 is not in use: `lsof -i :8000`

### Frontend can't connect to backend
- Check CORS settings in `backend/.env`
- Verify backend URL in `frontend/.env.local`
- Check firewall/network settings

### Pipeline fails
- Check Nextflow is installed: `nextflow -version`
- Verify all reference files exist
- Check Docker is running: `docker ps`
- Review pipeline logs in results directory

### File upload fails
- Check file size limits (adjust in backend if needed)
- Verify upload directory exists and is writable
- Check disk space: `df -h`

---

## File Structure

```
WholeExome/
├── backend/
│   ├── main.py              # FastAPI application
│   ├── auth.py              # Authentication logic
│   ├── database.py          # SQLAlchemy models
│   ├── schemas.py           # Pydantic schemas
│   ├── pipeline.py          # Nextflow pipeline runner
│   ├── config.py            # Configuration
│   ├── requirements.txt     # Python dependencies
│   └── .env                 # Environment variables
├── frontend/
│   ├── app/
│   │   ├── page.tsx         # Main page
│   │   ├── layout.tsx       # Layout wrapper
│   │   └── globals.css      # Global styles
│   ├── components/
│   │   ├── LoginForm.tsx    # Login component
│   │   ├── RegisterForm.tsx # Registration component
│   │   ├── Dashboard.tsx    # Main dashboard
│   │   ├── UploadForm.tsx   # File upload form
│   │   └── JobList.tsx      # Job list & downloads
│   ├── lib/
│   │   └── api.ts           # API client
│   └── package.json         # Node dependencies
├── uploads/                 # User uploaded files
├── results/                 # Pipeline results
└── main.nf                  # Nextflow pipeline
```

---

## Maintenance

### Database Backup
```bash
# Backup SQLite database
cp backend/wes_pipeline.db backend/wes_pipeline.db.backup
```

### Clean Old Results
```bash
# Remove results older than 30 days
find results/ -type f -mtime +30 -delete
find uploads/ -type f -mtime +30 -delete
```

### Monitor Disk Space
```bash
# Check space usage
du -sh uploads/ results/
df -h
```

---

## Support

For issues or questions:
- Check API logs: `backend/` directory
- Review Nextflow logs: `results/*/logs/`
- Original pipeline author: Robin Tomar (Aiimsgenomics@gmail.com)

---

## License

This pipeline maintains the original author attribution as required.
Author: Robin Tomar
GitHub: https://github.com/imrobintomar/Whole-Exome-Sequencing-Nextflow.git
