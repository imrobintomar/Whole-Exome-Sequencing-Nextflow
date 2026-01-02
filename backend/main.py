from fastapi import FastAPI, Depends, HTTPException, status, UploadFile, File, Form, BackgroundTasks
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from sqlalchemy.orm import Session
from pathlib import Path
from contextlib import asynccontextmanager
import uuid
import shutil
import os

from database import get_db, init_db, User, Job, JobStatus
from firebase_auth import get_current_user
from schemas import JobResponse, JobDetailResponse
from config import settings
from pipeline import run_job_async

# Lifespan event handler
@asynccontextmanager
async def lifespan(app: FastAPI):
    # Startup
    init_db()
    Path(settings.UPLOAD_DIR).mkdir(parents=True, exist_ok=True)
    Path(settings.RESULTS_DIR).mkdir(parents=True, exist_ok=True)
    yield
    # Shutdown (cleanup if needed)

app = FastAPI(title="WES Pipeline API", version="1.0.0", lifespan=lifespan)

# CORS configuration
app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.get_cors_origins(),
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Health check endpoint
@app.get("/")
def read_root():
    return {"status": "ok", "message": "WES Pipeline API is running"}

# Submit new job with FASTQ files
@app.post("/jobs/submit", response_model=JobResponse)
async def submit_job(
    background_tasks: BackgroundTasks,
    sample_name: str = Form(...),
    fastq_r1: UploadFile = File(...),
    fastq_r2: UploadFile = File(...),
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    # Validate file extensions
    if not fastq_r1.filename.endswith(('.fastq.gz', '.fq.gz')):
        raise HTTPException(status_code=400, detail="R1 file must be FASTQ (.fastq.gz or .fq.gz)")
    if not fastq_r2.filename.endswith(('.fastq.gz', '.fq.gz')):
        raise HTTPException(status_code=400, detail="R2 file must be FASTQ (.fastq.gz or .fq.gz)")

    # Generate unique job ID
    job_id = str(uuid.uuid4())

    # Create upload directory for this job
    upload_dir = Path(settings.UPLOAD_DIR) / job_id
    upload_dir.mkdir(parents=True, exist_ok=True)

    # Save uploaded files
    r1_path = upload_dir / f"{sample_name}_R1.fastq.gz"
    r2_path = upload_dir / f"{sample_name}_R2.fastq.gz"

    with open(r1_path, "wb") as f:
        shutil.copyfileobj(fastq_r1.file, f)
    with open(r2_path, "wb") as f:
        shutil.copyfileobj(fastq_r2.file, f)

    # Create job record in database
    job = Job(
        job_id=job_id,
        user_id=current_user.id,
        sample_name=sample_name,
        status=JobStatus.PENDING,
        fastq_r1_path=str(r1_path),
        fastq_r2_path=str(r2_path)
    )
    db.add(job)
    db.commit()
    db.refresh(job)

    # Start pipeline in background
    background_tasks.add_task(run_job_async, job_id, db)

    return job

# Get all jobs for current user
@app.get("/jobs", response_model=list[JobResponse])
def get_user_jobs(
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    jobs = db.query(Job).filter(Job.user_id == current_user.id).order_by(Job.created_at.desc()).all()
    return jobs

# Get job details
@app.get("/jobs/{job_id}", response_model=JobDetailResponse)
def get_job(
    job_id: str,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    job = db.query(Job).filter(Job.job_id == job_id, Job.user_id == current_user.id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    return job

# Download result file
@app.get("/jobs/{job_id}/download/{file_type}")
def download_file(
    job_id: str,
    file_type: str,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    job = db.query(Job).filter(Job.job_id == job_id, Job.user_id == current_user.id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    if job.status != JobStatus.COMPLETED:
        raise HTTPException(status_code=400, detail="Job not completed yet")

    # Map file type to path
    file_map = {
        "bam": job.bam_path,
        "raw_vcf": job.raw_vcf_path,
        "annotated_vcf": job.annotated_vcf_path,
        "filtered_tsv": job.filtered_tsv_path
    }

    file_path = file_map.get(file_type)
    if not file_path or not os.path.exists(file_path):
        raise HTTPException(status_code=404, detail=f"{file_type} file not found")

    return FileResponse(
        path=file_path,
        filename=os.path.basename(file_path),
        media_type="application/octet-stream"
    )

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
