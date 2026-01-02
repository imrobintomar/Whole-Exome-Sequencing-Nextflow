# Quick Start Guide

Get your WES Pipeline with web frontend up and running in minutes!

## Prerequisites

- Python 3.9+
- Node.js 18+
- Nextflow installed
- Docker installed
- All reference genomes and databases in place

---

## Backend Setup (5 minutes)

### 1. Configure Environment

```bash
cd backend
cp .env.example .env
nano .env  # Update paths to match your system
```

### 2. Start Backend

```bash
./run.sh
```

The backend will:
- Create virtual environment
- Install dependencies
- Start API server on port 8000

### 3. Test Backend

Open browser: `http://localhost:8000`

You should see: `{"status":"ok","message":"WES Pipeline API is running"}`

API Docs: `http://localhost:8000/docs`

---

## Frontend Setup (5 minutes)

### Option A: Test Locally First

```bash
cd frontend
npm install
cp .env.local.example .env.local
nano .env.local  # Set NEXT_PUBLIC_API_URL=http://localhost:8000
npm run dev
```

Visit: `http://localhost:3000`

### Option B: Deploy to Vercel Directly

1. **Expose your backend to internet:**

   ```bash
   # Install ngrok: https://ngrok.com/download
   ngrok http 8000
   ```

   Copy the HTTPS URL (e.g., `https://abc123.ngrok.io`)

2. **Update backend CORS:**

   Edit `backend/.env`:
   ```bash
   CORS_ORIGINS=http://localhost:3000,https://abc123.ngrok.io,https://your-app.vercel.app
   ```

   Restart backend: `./run.sh`

3. **Deploy to Vercel:**

   ```bash
   # Install Vercel CLI
   npm i -g vercel

   cd frontend
   vercel
   ```

   Follow prompts:
   - Set root directory: `frontend`
   - Add environment variable:
     - `NEXT_PUBLIC_API_URL` = `https://abc123.ngrok.io`

4. **Test it:**

   Visit your Vercel URL and create an account!

---

## First Usage

1. **Register Account**
   - Open your frontend URL
   - Click "Register here"
   - Create account

2. **Upload Files**
   - Login with your credentials
   - Go to "Upload Files" tab
   - Enter sample name
   - Upload R1 and R2 FASTQ files
   - Click "Submit Analysis"

3. **Monitor Progress**
   - Switch to "My Jobs" tab
   - Watch status change: PENDING → RUNNING → COMPLETED
   - Auto-refreshes every 10 seconds

4. **Download Results**
   - When COMPLETED, click download buttons:
     - BAM file
     - Raw VCF
     - Annotated VCF
     - Filtered TSV

---

## File Locations

After running the backend, files will be organized as:

```
WholeExome/
├── uploads/
│   └── {job-id}/
│       ├── {sample}_R1.fastq.gz
│       └── {sample}_R2.fastq.gz
├── results/
│   └── {job-id}/
│       ├── input/
│       └── output/
│           ├── Mapsam/
│           │   └── *_recal.bam
│           └── Germline_VCF/
│               ├── *.vcf.gz (raw)
│               ├── *_annovar_annotated.vcf
│               └── *_filtered.tsv
└── backend/
    └── wes_pipeline.db  # SQLite database
```

---

## Common Commands

### Backend

```bash
# Start backend
cd backend && ./run.sh

# View API docs
open http://localhost:8000/docs

# Check database
sqlite3 wes_pipeline.db "SELECT * FROM jobs;"

# View logs
tail -f nohup.out  # If running in background
```

### Frontend

```bash
# Local development
cd frontend && npm run dev

# Build for production
npm run build

# Deploy to Vercel
vercel --prod
```

### Maintenance

```bash
# Check disk usage
du -sh uploads/ results/

# Clean old files (30+ days)
find uploads/ -mtime +30 -type f -delete
find results/ -mtime +30 -type f -delete

# Backup database
cp backend/wes_pipeline.db backend/backup_$(date +%Y%m%d).db
```

---

## Troubleshooting

### "Port 8000 already in use"

```bash
# Find and kill process
lsof -ti:8000 | xargs kill -9

# Or use different port
# Edit backend/main.py, change port in uvicorn.run()
```

### "Frontend can't connect to backend"

```bash
# Check backend is running
curl http://localhost:8000

# Check CORS settings in backend/.env
# Check API URL in frontend/.env.local

# Verify ngrok is running (if using)
curl https://your-ngrok-url.ngrok.io
```

### "Pipeline not running"

```bash
# Check Nextflow installed
nextflow -version

# Check Docker running
docker ps

# Check reference files exist
ls -lh /media/drprabudh/m1/hg38/hg38.fa

# View Nextflow logs
cat results/{job-id}/output/logs/*.log
```

### "File upload fails"

```bash
# Check disk space
df -h

# Check upload directory permissions
ls -ld uploads/

# Make writable if needed
chmod 755 uploads/
```

---

## Production Tips

1. **Use a proper domain** instead of ngrok for production
   - Set up reverse proxy with Nginx
   - Use Let's Encrypt for SSL

2. **Run backend as systemd service** for auto-restart
   ```bash
   # Create /etc/systemd/system/wes-api.service
   sudo systemctl enable wes-api
   sudo systemctl start wes-api
   ```

3. **Set up monitoring**
   - Monitor disk space
   - Monitor job queue
   - Alert on failures

4. **Regular backups**
   - Backup database daily
   - Archive old results
   - Keep critical data redundant

5. **Security hardening**
   - Change default SECRET_KEY
   - Use strong passwords
   - Rate limit API endpoints
   - Regular security updates

---

## Next Steps

- Read full documentation: [README_FRONTEND.md](./README_FRONTEND.md)
- Customize frontend UI in `frontend/components/`
- Add email notifications in `backend/pipeline.py`
- Set up monitoring dashboard
- Configure automated backups

---

## Support

- Backend API docs: `http://localhost:8000/docs`
- Original pipeline: https://github.com/imrobintomar/Whole-Exome-Sequencing-Nextflow.git
- Check logs in `results/*/logs/` for pipeline issues
