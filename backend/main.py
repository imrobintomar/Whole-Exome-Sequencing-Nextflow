from dotenv import load_dotenv

# Load environment variables from .env file FIRST
load_dotenv()

from fastapi import FastAPI, Depends, HTTPException, status, UploadFile, File, Form, BackgroundTasks, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse
from sqlalchemy.orm import Session
from pathlib import Path
from contextlib import asynccontextmanager
import uuid
import shutil
import os
import traceback

from database import get_db, init_db, User, Job, JobStatus
from firebase_auth import get_current_user
from schemas import JobResponse, JobDetailResponse, GeneListFilter
from config import settings
from pipeline import run_job_async, PipelineRunner
from gene_panels import GenePanelManager
from acmg_classifier import ACMGClassifier
from constraint_data import get_constraint_db
from variant_analyzer import VariantAnalyzer
import threading
import pandas as pd

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

# Import SaaS modules (if available)
try:
    from modules.admin.routes import router as admin_router
    from modules.billing.routes import router as billing_router
    from modules.chat.routes import router as chat_router

    app.include_router(admin_router)
    app.include_router(billing_router)
    app.include_router(chat_router)
    print("✅ SaaS modules loaded successfully")
except Exception as e:
    print(f"⚠️  SaaS modules not available: {type(e).__name__}: {e}")
    import traceback
    traceback.print_exc()

# Middleware to handle ngrok-skip-browser-warning
@app.middleware("http")
async def add_ngrok_headers(request, call_next):
    response = await call_next(request)
    # Add headers for ngrok compatibility
    response.headers["ngrok-skip-browser-warning"] = "true"
    return response

# CORS configuration - must be added AFTER custom middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.get_cors_origins(),
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
    expose_headers=["*"],
)

# Global exception handler to ensure CORS headers are sent even on errors
def get_analysis_file(job: Job) -> str:
    """
    Get the analysis file path with fallback logic.

    Priority:
    1. filtered_tsv_path (Final_.txt with UniqueID) - new pipeline
    2. annotated_vcf_path (ANNOVAR multianno.txt) - old pipeline or if Final_.txt missing

    For old jobs, automatically add UniqueID column if missing.
    """
    # Try filtered_tsv_path first (new pipeline output)
    if job.filtered_tsv_path and os.path.exists(job.filtered_tsv_path):
        return job.filtered_tsv_path

    # Fallback to ANNOVAR file (old pipeline or incomplete pipeline)
    if job.annotated_vcf_path and os.path.exists(job.annotated_vcf_path):
        # Check if this ANNOVAR file already has UniqueID column
        try:
            df = pd.read_csv(job.annotated_vcf_path, sep='\t', nrows=1, encoding='utf-8', on_bad_lines='skip')
            if 'UniqueID' in df.columns:
                # Already has UniqueID, can use directly
                return job.annotated_vcf_path

            # Need to add UniqueID - create temp file with UniqueID added
            print(f"⚠️  Job {job.job_id}: Using ANNOVAR file without UniqueID, adding it now...")
            return job.annotated_vcf_path  # For now, return as-is and handle in VariantAnalyzer
        except Exception as e:
            print(f"⚠️  Error checking UniqueID column: {e}")
            return job.annotated_vcf_path

    return None

@app.exception_handler(Exception)
async def global_exception_handler(request: Request, exc: Exception):
    """
    Catch all exceptions and return proper CORS-enabled error responses
    This ensures CORS headers are added even when endpoints fail
    """
    # Log the full exception for debugging
    print(f"❌ Unhandled exception in {request.method} {request.url.path}:")
    print(f"   Error: {type(exc).__name__}: {str(exc)}")
    traceback.print_exc()

    # Determine status code
    if isinstance(exc, HTTPException):
        status_code = exc.status_code
        detail = exc.detail
    else:
        status_code = 500
        detail = f"Internal server error: {str(exc)}"

    # Return JSON response with CORS headers
    return JSONResponse(
        status_code=status_code,
        content={"detail": detail},
        headers={
            "Access-Control-Allow-Origin": request.headers.get("origin", "*"),
            "Access-Control-Allow-Credentials": "true",
            "Access-Control-Allow-Methods": "*",
            "Access-Control-Allow-Headers": "*",
        }
    )

# Health check endpoint
@app.get("/")
def read_root():
    return {"status": "ok", "message": "WES Pipeline API is running"}

# Submit new job WITH billing/usage enforcement (SaaS wrapper)
@app.post("/jobs/submit-with-billing", response_model=JobResponse)
async def submit_job_with_billing(
    background_tasks: BackgroundTasks,
    sample_name: str = Form(...),
    fastq_r1: UploadFile = File(...),
    fastq_r2: UploadFile = File(...),
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Submit job with usage limit enforcement (wrapper for SaaS)
    This is the recommended endpoint for production use
    """
    try:
        from middleware.subscription_guard import enforce_usage_limit
        from services.usage_tracking_service import UsageTrackingService
        from services.audit_service import AuditService

        # Enforce usage limits
        try:
            billing_info = await enforce_usage_limit(current_user)
        except HTTPException as e:
            # Usage limit exceeded or payment required
            raise e

        # Proceed with job submission (same logic as original endpoint)
        job = await submit_job(
            background_tasks=background_tasks,
            sample_name=sample_name,
            fastq_r1=fastq_r1,
            fastq_r2=fastq_r2,
            current_user=current_user,
            db=db
        )

        # Increment usage counter
        usage_service = UsageTrackingService(db)
        usage_service.increment_job_count(current_user.uid)

        # Audit log
        audit = AuditService(db)
        audit.log_action(
            action="job_submitted",
            user_id=current_user.uid,
            resource_type="job",
            resource_id=job.job_id,
            metadata={
                "sample_name": sample_name,
                "plan": billing_info["plan"].name if billing_info.get("plan") else "Free"
            }
        )

        return job

    except ImportError:
        # SaaS modules not available, fall back to regular submission
        return await submit_job(
            background_tasks=background_tasks,
            sample_name=sample_name,
            fastq_r1=fastq_r1,
            fastq_r2=fastq_r2,
            current_user=current_user,
            db=db
        )


# Submit new job with FASTQ files (original endpoint - kept for backward compatibility)
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

    # Start pipeline in background (don't pass db session)
    background_tasks.add_task(run_job_async, job_id)

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

# Pipeline Control Endpoints
@app.post("/jobs/{job_id}/cancel")
def cancel_job(
    job_id: str,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """Cancel a running pipeline job"""
    job = db.query(Job).filter(Job.job_id == job_id, Job.user_id == current_user.id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    try:
        runner = PipelineRunner(job, db)
        success = runner.cancel_pipeline()
        return {"message": "Pipeline cancelled successfully" if success else "Pipeline already completed", "job_id": job_id}
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/jobs/{job_id}/rerun")
def rerun_job(
    job_id: str,
    background_tasks: BackgroundTasks,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """Rerun a pipeline from scratch"""
    job = db.query(Job).filter(Job.job_id == job_id, Job.user_id == current_user.id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    # Run pipeline in background thread with proper session management
    def run_rerun():
        from database import SessionLocal
        db_thread = SessionLocal()
        try:
            # Re-fetch job in this thread's session
            job_thread = db_thread.query(Job).filter(Job.job_id == job_id).first()
            if job_thread:
                runner = PipelineRunner(job_thread, db_thread)
                runner.rerun_pipeline()
        except Exception as e:
            print(f"Error in rerun thread: {e}")
        finally:
            db_thread.close()

    thread = threading.Thread(target=run_rerun)
    thread.daemon = True
    thread.start()

    return {"message": "Pipeline rerun started", "job_id": job_id}

@app.post("/jobs/{job_id}/resume")
def resume_job(
    job_id: str,
    background_tasks: BackgroundTasks,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """Resume a failed pipeline"""
    job = db.query(Job).filter(Job.job_id == job_id, Job.user_id == current_user.id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    if job.status != JobStatus.FAILED:
        raise HTTPException(status_code=400, detail=f"Cannot resume job with status: {job.status}")

    # Run pipeline in background thread with proper session management
    def run_resume():
        from database import SessionLocal
        db_thread = SessionLocal()
        try:
            # Re-fetch job in this thread's session
            job_thread = db_thread.query(Job).filter(Job.job_id == job_id).first()
            if job_thread:
                runner = PipelineRunner(job_thread, db_thread)
                runner.resume_pipeline()
        except Exception as e:
            print(f"Error in resume thread: {e}")
        finally:
            db_thread.close()

    thread = threading.Thread(target=run_resume)
    thread.daemon = True
    thread.start()

    return {"message": "Pipeline resume started", "job_id": job_id}

@app.delete("/jobs/{job_id}")
def delete_job(
    job_id: str,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """Delete a job and its associated files (user can only delete their own jobs)"""
    job = db.query(Job).filter(Job.job_id == job_id, Job.user_id == current_user.id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    # Cancel if running
    if job.status == JobStatus.RUNNING:
        try:
            runner = PipelineRunner(job, db)
            runner.cancel_pipeline()
        except Exception as e:
            print(f"Warning: Failed to cancel running job before deletion: {e}")

    # Delete associated files
    try:
        # Delete upload directory
        upload_dir = Path(settings.UPLOAD_DIR) / job_id
        if upload_dir.exists():
            shutil.rmtree(upload_dir)

        # Delete results directory
        results_dir = Path(settings.RESULTS_DIR) / job_id
        if results_dir.exists():
            shutil.rmtree(results_dir)
    except Exception as e:
        print(f"Warning: Failed to delete some files for job {job_id}: {e}")

    # Delete database record
    db.delete(job)
    db.commit()

    return {"message": "Job deleted successfully", "job_id": job_id}

# Check if BAM file exists
@app.get("/jobs/{job_id}/files/bam")
def check_bam_file(
    job_id: str,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """Check if BAM file exists for this job"""
    job = db.query(Job).filter(Job.job_id == job_id, Job.user_id == current_user.id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    if job.status != JobStatus.COMPLETED:
        raise HTTPException(status_code=404, detail="Job not completed yet")

    if not job.bam_path or not os.path.exists(job.bam_path):
        raise HTTPException(status_code=404, detail="BAM file not found")

    return {"available": True, "path": job.bam_path}

# Check if VCF file exists
@app.get("/jobs/{job_id}/files/vcf")
def check_vcf_file(
    job_id: str,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """Check if VCF file exists for this job"""
    job = db.query(Job).filter(Job.job_id == job_id, Job.user_id == current_user.id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    if job.status != JobStatus.COMPLETED:
        raise HTTPException(status_code=404, detail="Job not completed yet")

    if not job.raw_vcf_path or not os.path.exists(job.raw_vcf_path):
        raise HTTPException(status_code=404, detail="VCF file not found")

    return {"available": True, "path": job.raw_vcf_path}

# Download result file
@app.get("/jobs/{job_id}/download/{file_type}")
def download_file(
    job_id: str,
    file_type: str,
    index: bool = False,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    job = db.query(Job).filter(Job.job_id == job_id, Job.user_id == current_user.id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    if job.status != JobStatus.COMPLETED:
        raise HTTPException(status_code=400, detail="Job not completed yet")

    # Map file type to path and custom filename (updated for simplified pipeline)
    file_map = {
        "bam": (job.bam_path, f"{job.sample_name}_recall.bam"),
        "bam.bai": (f"{job.bam_path}.bai" if job.bam_path else None, f"{job.sample_name}_recall.bam.bai"),
        "raw_vcf": (job.raw_vcf_path, f"{job.sample_name}.vcf.gz"),
        "raw_vcf.tbi": (f"{job.raw_vcf_path}.tbi" if job.raw_vcf_path else None, f"{job.sample_name}.vcf.gz.tbi"),
        # annotated_vcf now points to ANNOVAR TXT file
        "annotated_vcf": (job.annotated_vcf_path, f"{job.sample_name}.annovar.hg38_multianno.txt"),
        # filtered_tsv now points to Final_.txt (ANNOVAR format with UniqueID)
        "filtered_tsv": (job.filtered_tsv_path, f"{job.sample_name}_Final_.txt")
    }

    file_info = file_map.get(file_type)
    if not file_info:
        raise HTTPException(status_code=404, detail=f"Invalid file type: {file_type}")

    file_path, download_filename = file_info

    # If requesting index file via query parameter (backward compatibility)
    if index:
        if file_type == "bam":
            file_path = f"{file_path}.bai"
            download_filename = f"{download_filename}.bai"
        elif file_type == "raw_vcf":
            # Only raw_vcf has index, annotated_vcf is now TXT format
            file_path = f"{file_path}.tbi"
            download_filename = f"{download_filename}.tbi"
        else:
            raise HTTPException(status_code=400, detail=f"Index not available for {file_type}")

    if not file_path or not os.path.exists(file_path):
        raise HTTPException(status_code=404, detail=f"{file_type} {'index' if index else 'file'} not found")

    # Set appropriate media type
    media_type = "application/octet-stream"
    if file_type == "bam" or (index and file_type == "bam"):
        media_type = "application/octet-stream"
    elif file_type == "raw_vcf":
        media_type = "application/gzip" if not index else "application/octet-stream"
    elif file_type in ["annotated_vcf", "filtered_tsv"]:
        # Both are now text files
        media_type = "text/tab-separated-values"

    return FileResponse(
        path=file_path,
        filename=download_filename,
        media_type=media_type
    )

# Download filtered variants by gene list
@app.post("/jobs/{job_id}/download/filtered")
async def download_filtered_variants(
    job_id: str,
    gene_filter: GeneListFilter,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """Download TSV file filtered to specific genes (validated gene list)"""
    import pandas as pd
    import tempfile

    job = db.query(Job).filter(Job.job_id == job_id, Job.user_id == current_user.id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    if job.status != JobStatus.COMPLETED:
        raise HTTPException(status_code=400, detail="Job not completed yet")

    # Get analysis file with fallback to ANNOVAR
    analysis_file = get_analysis_file(job)
    if not analysis_file:
        raise HTTPException(status_code=404, detail="TSV file not found")

    # Extract validated genes from request
    genes = gene_filter.genes

    # Read and filter TSV
    try:
        # Handle encoding issues from pipeline output
        try:
            df = pd.read_csv(analysis_file, sep='\t', encoding='utf-8')
        except UnicodeDecodeError:
            df = pd.read_csv(analysis_file, sep='\t', encoding='latin1', on_bad_lines='skip')

        # Filter by genes (handle multiple genes per row)
        def gene_in_list(gene_str):
            if pd.isna(gene_str):
                return False
            variant_genes = [g.strip().upper() for g in str(gene_str).split(',')]
            return any(g in genes for g in variant_genes)

        gene_column = 'Gene.refGeneWithVer'
        if gene_column in df.columns:
            filtered_df = df[df[gene_column].apply(gene_in_list)]
        else:
            raise HTTPException(status_code=500, detail="Gene column not found in TSV")

        # Write to temporary file with proper cleanup
        temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False)
        try:
            filtered_df.to_csv(temp_file.name, sep='\t', index=False)
            temp_path = temp_file.name
            temp_file.close()

            # Safe cleanup function
            def cleanup_temp_file():
                try:
                    if os.path.exists(temp_path):
                        os.unlink(temp_path)
                except Exception as e:
                    print(f"Warning: Failed to cleanup temp file {temp_path}: {e}")

            return FileResponse(
                path=temp_path,
                filename=f"{job.sample_name}_filtered_{len(genes)}genes.tsv",
                media_type="text/tab-separated-values",
                background=cleanup_temp_file
            )
        except Exception as write_error:
            # Clean up temp file if writing failed
            temp_file.close()
            if os.path.exists(temp_file.name):
                os.unlink(temp_file.name)
            raise HTTPException(status_code=500, detail=f"Failed to write filtered variants: {str(write_error)}")

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to filter variants: {str(e)}")

# Gene Panel Endpoints
@app.get("/panels/search")
def search_gene_panels(query: str, current_user: User = Depends(get_current_user)):
    """Search for gene panels from PanelApp"""
    panels = GenePanelManager.search_panels(query)
    return {"results": panels}

@app.get("/panels/{panel_id}/genes")
def get_panel_genes(
    panel_id: int,
    confidence_level: str = "3",
    current_user: User = Depends(get_current_user)
):
    """Get genes from a specific panel"""
    genes = GenePanelManager.get_panel_genes(panel_id, confidence_level)
    return {"panel_id": panel_id, "genes": genes, "count": len(genes)}

@app.get("/panels/acmg-sf")
def get_acmg_secondary_findings(current_user: User = Depends(get_current_user)):
    """Get ACMG Secondary Findings v3.2 gene list"""
    genes = GenePanelManager.get_acmg_secondary_findings_genes()
    return {"genes": genes, "count": len(genes), "version": "3.2"}

@app.post("/jobs/{job_id}/apply-panel")
async def apply_gene_panel(
    job_id: str,
    gene_filter: GeneListFilter,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Apply gene panel filter to job's TSV and return filtered variants with statistics
    """
    import pandas as pd

    job = db.query(Job).filter(Job.job_id == job_id, Job.user_id == current_user.id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    if job.status != JobStatus.COMPLETED:
        raise HTTPException(status_code=400, detail="Job not completed yet")

    # Get analysis file with fallback to ANNOVAR
    analysis_file = get_analysis_file(job)
    if not analysis_file:
        raise HTTPException(status_code=404, detail="TSV file not found")

    try:
        # Read TSV (handle encoding issues from pipeline output)
        try:
            df = pd.read_csv(analysis_file, sep='\t', encoding='utf-8')
        except UnicodeDecodeError:
            df = pd.read_csv(analysis_file, sep='\t', encoding='latin1', on_bad_lines='skip')

        genes = gene_filter.genes

        # Filter by genes using GenePanelManager
        filtered_df = GenePanelManager.filter_variants_by_genes(df, genes)

        # Calculate statistics
        total_variants = len(filtered_df)

        # Chromosome distribution
        chr_dist = {}
        if 'Chr' in filtered_df.columns:
            chr_counts = filtered_df['Chr'].value_counts().to_dict()
            chr_dist = {str(k): int(v) for k, v in chr_counts.items()}

        # Gene distribution (top 20)
        gene_dist = {}
        gene_column = 'Gene.refGeneWithVer'
        if gene_column in filtered_df.columns:
            gene_counts = filtered_df[gene_column].value_counts().head(20).to_dict()
            gene_dist = {str(k): int(v) for k, v in gene_counts.items()}

        # Functional impact
        func_impact = {}
        if 'Func.refGeneWithVer' in filtered_df.columns:
            func_counts = filtered_df['Func.refGeneWithVer'].value_counts().to_dict()
            func_impact = {str(k): int(v) for k, v in func_counts.items()}

        # Exonic subcategories
        exonic_subcat = {}
        if 'ExonicFunc.refGeneWithVer' in filtered_df.columns:
            exonic_counts = filtered_df['ExonicFunc.refGeneWithVer'].value_counts().to_dict()
            exonic_subcat = {str(k): int(v) for k, v in exonic_counts.items()}

        # Clinical significance
        clinical_sig = {}
        if 'CLNSIG' in filtered_df.columns:
            clinical_counts = filtered_df['CLNSIG'].value_counts().to_dict()
            clinical_sig = {str(k): int(v) for k, v in clinical_counts.items()}

        # Convert filtered variants to list of dicts (limit to first 1000 for performance)
        variants_list = filtered_df.head(1000).to_dict('records')

        return {
            "job_id": job_id,
            "sample_name": job.sample_name,
            "applied_genes": genes,
            "gene_count": len(genes),
            "total_variants": total_variants,
            "showing_variants": min(1000, total_variants),
            "variants": variants_list,
            "statistics": {
                "chromosome_distribution": chr_dist,
                "gene_distribution": gene_dist,
                "functional_impact": func_impact,
                "exonic_subcategories": exonic_subcat,
                "clinical_significance": clinical_sig
            }
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to apply gene panel: {str(e)}")

# ACMG Classification Endpoint
@app.post("/classify/acmg")
async def classify_variant_acmg(
    variant: dict,
    current_user: User = Depends(get_current_user)
):
    """
    Classify a variant using ACMG/AMP 2015 guidelines

    Request body should contain variant information:
    {
        "consequence": "missense_variant",
        "gene": "BRCA1",
        "af_gnomad": 0.00001,
        "cadd_phred": 25,
        "revel_score": 0.7,
        "sift_pred": "deleterious",
        "polyphen_pred": "probably_damaging",
        "clinvar_sig": "Pathogenic",
        "spliceai_max": 0.1
    }
    """
    try:
        # Add gene constraint data
        gene = variant.get("gene", "")
        if gene:
            constraint_db = get_constraint_db()
            constraint = constraint_db.get_gene_constraint(gene)
            variant["pli"] = constraint.get("pli", 0)
            variant["loeuf"] = constraint.get("loeuf", 1.0)

        # Classify variant
        classifier = ACMGClassifier()
        result = classifier.classify_variant(variant)

        return result

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Classification failed: {str(e)}")

@app.post("/jobs/{job_id}/classify")
async def classify_job_variants(
    job_id: str,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Apply ACMG classification to all variants in a job's TSV file
    """
    import pandas as pd

    job = db.query(Job).filter(Job.job_id == job_id, Job.user_id == current_user.id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    if job.status != JobStatus.COMPLETED:
        raise HTTPException(status_code=400, detail="Job not completed yet")

    # Get analysis file with fallback to ANNOVAR
    analysis_file = get_analysis_file(job)
    if not analysis_file:
        raise HTTPException(status_code=404, detail="TSV file not found")

    try:
        # Read TSV (handle encoding issues from pipeline output)
        try:
            df = pd.read_csv(analysis_file, sep='\t', encoding='utf-8')
        except UnicodeDecodeError:
            # Fallback to latin1 encoding if UTF-8 fails
            df = pd.read_csv(analysis_file, sep='\t', encoding='latin1', on_bad_lines='skip')

        # Get constraint DB
        constraint_db = get_constraint_db()

        # Classify each variant
        classifications = []
        classifier = ACMGClassifier()

        for idx, row in df.iterrows():
            # Helper function to safely get values
            def safe_get(column_name, default=""):
                if column_name in df.columns:
                    val = row[column_name]
                    if pd.isna(val):
                        return default
                    return val
                return default

            def safe_float(column_name, default=0.0):
                if column_name in df.columns:
                    val = row[column_name]
                    if pd.isna(val) or val == "" or val == ".":
                        return default
                    try:
                        return float(val)
                    except (ValueError, TypeError):
                        return default
                return default

            # Map TSV columns to classifier input
            variant = {
                "consequence": safe_get("ExonicFunc.refGeneWithVer", ""),
                "gene": safe_get("Gene.refGeneWithVer", ""),
                "af_gnomad": safe_float("AF", 0.0),
                "cadd_phred": safe_float("CADD_phred", 0.0),
                "revel_score": safe_float("REVEL_score", 0.0),
                "sift_pred": safe_get("SIFT_pred", ""),
                "polyphen_pred": safe_get("Polyphen2_HDIV_pred", ""),
                "clinvar_sig": safe_get("CLNSIG", ""),  # May not exist in TSV
            }

            # Add constraint data
            gene = variant["gene"]
            if gene and gene != ".":
                constraint = constraint_db.get_gene_constraint(gene)
                variant["pli"] = constraint.get("pli", 0)
                variant["loeuf"] = constraint.get("loeuf", 1.0)

            # Classify
            result = classifier.classify_variant(variant)

            # Get UniqueID with fallback logic
            # Priority:
            # 1. UniqueID column (from Final_.txt with add_unique_id.py)
            # 2. Construct from ANNOVAR columns: Chr:Start:Ref:Alt
            # 3. Construct from VCF columns: CHROM:POS:REF:ALT (legacy)
            unique_id = safe_get('UniqueID', '')
            if not unique_id or unique_id == '.':
                # Try ANNOVAR format
                chr_val = safe_get('Chr', safe_get('CHROM', ''))
                start_val = safe_get('Start', safe_get('POS', ''))
                ref_val = safe_get('Ref', safe_get('REF', ''))
                alt_val = safe_get('Alt', safe_get('ALT', ''))
                unique_id = f"{chr_val}:{start_val}:{ref_val}:{alt_val}"

            classifications.append({
                "position": unique_id,
                "gene": gene,
                "consequence": variant["consequence"],
                "classification": result["classification"],
                "evidence_summary": result["evidence_summary"],
                "met_criteria": result["met_criteria"]
            })

        return {
            "job_id": job_id,
            "sample_name": job.sample_name,
            "total_variants": len(classifications),
            "classifications": classifications,
            "summary": {
                "pathogenic": sum(1 for c in classifications if c["classification"] == "Pathogenic"),
                "likely_pathogenic": sum(1 for c in classifications if c["classification"] == "Likely Pathogenic"),
                "vus": sum(1 for c in classifications if c["classification"] == "Uncertain Significance"),
                "likely_benign": sum(1 for c in classifications if c["classification"] == "Likely Benign"),
                "benign": sum(1 for c in classifications if c["classification"] == "Benign"),
            }
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Classification failed: {str(e)}")

# Variant Visualization Endpoint
@app.get("/jobs/{job_id}/variant-metrics")
async def get_variant_metrics(
    job_id: str,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Get variant analysis metrics for visualization
    """
    job = db.query(Job).filter(Job.job_id == job_id, Job.user_id == current_user.id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    if job.status != JobStatus.COMPLETED:
        raise HTTPException(status_code=400, detail="Job not completed yet")

    # Get analysis file with fallback to ANNOVAR
    analysis_file = get_analysis_file(job)
    if not analysis_file:
        raise HTTPException(status_code=404, detail="TSV file not found")

    try:
        analyzer = VariantAnalyzer(analysis_file)
        metrics = analyzer.get_all_metrics()

        return {
            "job_id": job_id,
            "sample_name": job.sample_name,
            "metrics": metrics
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Analysis failed: {str(e)}")

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
